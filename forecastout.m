% Clear workspace and command window
clear;
clc;

% Define file path
filledDataPath = 'E:/OneDrive - Lancaster University/Msc MBF/2. ECO420 Dissertation/Week 3 Data/Final data/Data for matlab/stationaryfilled.xlsx';

% Check if the file exists
if ~isfile(filledDataPath)
    error('File does not exist. Please check the file path.');
end

% Read filled data
dataMatrix = readmatrix(filledDataPath);
dataTable = readtable(filledDataPath, 'ReadVariableNames', true);

% Display the first few rows to confirm data reading
disp('First few rows of the data:');
disp(dataTable(1:5, :));

% Standardize data
dataMatrix = zscore(dataMatrix);

% Extract target variables
y_targets = dataMatrix(:, 1:3); % Assume the first 3 columns are target variables

% Define the number of years for in-sample and out-sample
in_sample_years = 10;
out_sample_years = 20;

% Define the number of observations per year (assumed to be monthly data)
observations_per_year = 12;

% Calculate the number of observations for in-sample and out-sample
in_sample_obs = in_sample_years * observations_per_year;
out_sample_obs = out_sample_years * observations_per_year;

% Check if the data matrix has enough observations
if size(dataMatrix, 1) < (in_sample_obs + out_sample_obs)
    error('Not enough observations in the data matrix for the specified in-sample and out-sample periods.');
end

% Split the data into in-sample and out-sample sets
data_in_sample = dataMatrix(1:in_sample_obs, :);
data_out_sample = dataMatrix((in_sample_obs + 1):(in_sample_obs + out_sample_obs), :);

% Extract target variables for in-sample and out-sample
y_in_sample = data_in_sample(:, 1:3);
y_out_sample = data_out_sample(:, 1:3);

% Extract predictor variables (assuming predictors are all columns except the first 3)
X_in_sample = data_in_sample(:, 4:end);
X_out_sample = data_out_sample(:, 4:end);

%%
% Perform PCA on the entire set of predictor variables (in-sample)
[coeff, score, latent, tsquared, explained, mu] = pca(X_in_sample, 'Rows', 'pairwise');

% Determine the number of principal components to use based on the cumulative explained variance
cumulative_explained = cumsum(explained);
num_components = find(cumulative_explained >= 95, 1); % Choose the number of components that explain at least 95% of the variance

% Use the selected principal components as latent factors
latent_factors_in_sample = score(:, 1:num_components);

%%
% Initialize result storage for predictions and evaluation metrics
predictions = zeros(out_sample_obs, size(y_targets, 2));
me_values = zeros(size(y_targets, 2), 1);
rmse_values = zeros(size(y_targets, 2), 1);
mae_values = zeros(size(y_targets, 2), 1);
mpe_values = zeros(size(y_targets, 2), 1);
mape_values = zeros(size(y_targets, 2), 1);
mase_values = zeros(size(y_targets, 2), 1);
acf1_values = zeros(size(y_targets, 2), 1);

% Initialize R-squared values storage
r_squared_values = zeros(size(y_targets, 2), 1);

% Perform regression analysis and predict out-sample data
for i = 1:size(y_targets, 2)
    y_it_in_sample = y_in_sample(:, i);

    % Perform ordinary least squares regression
    mdl = fitlm(latent_factors_in_sample, y_it_in_sample);

    % Extract R-squared value
    r_squared_values(i) = mdl.Rsquared.Ordinary;

    % Project the out-sample predictor variables onto the principal components
    latent_factors_out_sample = (X_out_sample - mu) * coeff(:, 1:num_components);

    % Predict the out-sample target variable
    predictions(:, i) = predict(mdl, latent_factors_out_sample);

    % Calculate residuals
    residuals = y_out_sample(:, i) - predictions(:, i);

    % Calculate evaluation metrics
    me_values(i) = mean(residuals);
    rmse_values(i) = sqrt(mean(residuals.^2));
    mae_values(i) = mean(abs(residuals));
    mpe_values(i) = mean((residuals ./ y_out_sample(:, i)) * 100);
    mape_values(i) = mean(abs(residuals ./ y_out_sample(:, i)) * 100);
    mase_values(i) = mae_values(i) / mean(abs(diff(y_in_sample(:, i))));

    % Calculate first lag autocorrelation
    acf_values = autocorr(residuals, 1);
    acf1_values(i) = acf_values(2); % The first element is the lag 0 autocorrelation which is always 1, so we take the second element
end

%%
% Create results table for the current group
r_squared_table = array2table(r_squared_values', 'VariableNames', {'SPGSEN Index', 'SPGSIN Index', 'SPGSAGS Index'});

% Display the results
disp('R-squared values:');
disp(r_squared_table);

%%
% Display evaluation metrics
disp('Evaluation Metrics:');
for i = 1:size(y_targets, 2)
    fprintf('Target %d:\n', i);
    fprintf('  ME    = %.4f\n', me_values(i));
    fprintf('  RMSE  = %.4f\n', rmse_values(i));
    fprintf('  MAE   = %.4f\n', mae_values(i));
    fprintf('  MPE   = %.4f\n', mpe_values(i));
    fprintf('  MAPE  = %.4f\n', mape_values(i));
    fprintf('  MASE  = %.4f\n', mase_values(i));
    fprintf('  ACF1  = %.4f\n', acf1_values(i));
end


% Define target names for better visualization
target_names = {'SPGSEN Index', 'SPGSIN Index', 'SPGSAGS Index'};

% Loop through each target variable and plot
for i = 1:size(y_targets, 2)
    figure;
    hold on;
    % Plot original data
    plot(y_out_sample(:, i), 'b', 'DisplayName', 'Original Data');
    % Plot predictions
    plot(predictions(:, i), 'r--', 'DisplayName', 'Predicted Data');
    hold off;
    
    % Set plot title and labels
    title(['Original vs Predicted Data for ', target_names{i}]);
    xlabel('Time');
    ylabel('Value');
    legend('show');
    
    % Save the figure
    saveas(gcf, ['Prediction_vs_Original_', target_names{i}, '.png']);
end
