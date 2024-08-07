% Clear workspace and close figures
clear;
clc;
close all;

% Load data from Excel files
monthlyTarget = readtable('MonthlyTarget.xlsx');
financialVars = readtable('Financial variables.xlsx');
transFRED = readtable('Trans FRED.xlsx');

% Remove the date columns (assuming the date column is the first one in each table)
monthlyTarget = monthlyTarget(:, 2:end);
financialVars = financialVars(:, 2:end);
transFRED = transFRED(:, 2:end);

% Combine tables (assuming all datasets now have the same number of rows in the same order)
combinedData = [monthlyTarget, financialVars, transFRED];

% Convert table to matrix
dataMatrix = table2array(combinedData); % All columns are now numerical data

% Check for NaNs and Infs
if any(isnan(dataMatrix(:))) || any(isinf(dataMatrix(:)))
    fprintf('Data contains NaNs or Infs. Proceeding with data cleaning...\n');
    
    % Replace NaNs with column mean (can also use interpolation or other methods)
    colMeans = nanmean(dataMatrix);
    dataMatrix(isnan(dataMatrix)) = repmat(colMeans(isnan(dataMatrix)), sum(isnan(dataMatrix(:))), 1);
    
    % Replace Infs with large finite values
    dataMatrix(isinf(dataMatrix)) = max(dataMatrix(~isinf(dataMatrix)), [], 'all') * 10;

    % Verify data after cleaning
    if any(isnan(dataMatrix(:))) || any(isinf(dataMatrix(:)))
        error('Data still contains NaNs or Infs after cleaning.');
    end
else
    fprintf('Data does not contain NaNs or Infs.\n');
end
%%
% Standardize the data
meanData = mean(dataMatrix);
stdData = std(dataMatrix);
stdData(stdData == 0) = 1; % Avoid division by zero
dataMatrix = (dataMatrix - meanData) ./ stdData;

% Check standardized data for NaNs and Infs
if any(isnan(dataMatrix(:))) || any(isinf(dataMatrix(:)))
    error('Standardized data contains NaNs or Infs.');
end

% Compute the cross-covariance matrix
T = size(dataMatrix, 1); % Number of time periods
N = size(dataMatrix, 2); % Number of variables
X_t = dataMatrix';
covMatrix = (X_t * X_t') / T;

% Check the covariance matrix for NaNs or Infs
if any(isnan(covMatrix(:))) || any(isinf(covMatrix(:)))
    error('Covariance matrix contains NaNs or Infs.');
end

% Perform eigenvalue decomposition
[eigVec, eigVal] = eig(covMatrix);
eigVal = diag(eigVal); % Extract diagonal elements of eigenvalue matrix
[eigVal, sortIdx] = sort(eigVal, 'descend');
eigVec = eigVec(:, sortIdx);

% Display eigenvalues
disp('Eigenvalues:');
disp(eigVal);

% Select the top 'r' eigenvectors
r = 3; % Number of factors (adjust as necessary)
if r > length(eigVal)
    error('Number of factors exceeds the number of available eigenvalues.');
end
V_r = eigVec(:, 1:r);
D_r = diag(eigVal(1:r));

% Estimate the factors
F_hat = V_r' * X_t;

% Check for NaNs in factors
if any(isnan(F_hat(:))) || any(isinf(F_hat(:)))
    error('Estimated factors contain NaNs or Infs.');
end

% Display first few estimated factors for inspection
disp('Estimated Factors (F_hat):');
disp(F_hat(1:5, :));
%%
% Compute the estimated covariance matrix of the common components
Gamma_chi_0 = V_r * D_r * V_r';

% Prepare for forecasting
% Assuming the target variable is the first column of dataMatrix
y_it = dataMatrix(:, 1);
X_t_lag = [y_it(1:end-1), X_t(:, 1:end-1)']; % Lagged values

% Fit the forecasting model
% Define the number of lags
numLags = 2;
y_lags = lagmatrix(y_it, 1:numLags);
y_lags = y_lags(numLags+1:end, :);
X_t_lag = X_t(numLags+1:end, :);

% Combine factors and lags for regression
regressionMatrix = [F_hat(:, numLags+1:end)', y_lags];
y = y_it(numLags+1:end);

% Perform Ordinary Least Squares regression
[b, bint, r, rint, stats] = regress(y, [ones(size(regressionMatrix, 1), 1), regressionMatrix]);

% Output regression results
fprintf('Regression Results:\n');
disp('Coefficients:');
disp(b);
disp('R-squared:');
disp(stats(1));

% Forecasting
h = 1; % Forecast horizon (adjust as necessary)

% Check dimensions of F_hat for the forecast horizon
if size(F_hat, 2) < numLags
    error('Insufficient number of factors for the forecast horizon.');
end

% Prepare the forecast matrix
% Ensure the number of columns in forecastMatrix matches the length of b
forecastMatrix = [ones(h, 1), F_hat(:, end-h+1:end)'];

% Print dimensions for debugging
disp('Forecast Matrix Dimensions:');
disp(size(forecastMatrix));
disp('Regression Coefficients Dimensions:');
disp(size(b));

% Ensure dimensions match
if size(forecastMatrix, 2) ~= length(b)
    error('Dimension mismatch: Forecast matrix and regression coefficients are not compatible.');
end

% Perform forecast calculation
forecast = forecastMatrix * b;

% Display forecast
fprintf('Forecast for next %d periods:\n', h);
disp(forecast);

% In-Sample and Out-of-Sample Analysis
% Here you would implement the specific analysis code based on your dataset and requirements.
% Example: Split data for out-of-sample analysis, calculate forecasting errors, etc.

% Example for out-of-sample analysis (pseudo-code):
% splitPoint = floor(0.8 * size(dataMatrix, 1));
% trainData = dataMatrix(1:splitPoint, :);
% testData = dataMatrix(splitPoint+1:end, :);
% 
% % Re-estimate model on trainData
% % Forecast on testData
% % Evaluate performance using metrics (e.g., RMSE)

