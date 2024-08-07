% Clear workspace and command window
clear;
clc;

% Define file path for the combined data file
combinedDataPath = 'E:/OneDrive - Lancaster University/Msc MBF/2. ECO420 Dissertation/Week 3 Data/Final data/Data for matlab/stationaryfilled.xlsx';

% Check if the file exists
if ~isfile(combinedDataPath)
    error('File does not exist. Please check the file path.');
end

% Read combined data
dataTable = readtable(combinedDataPath, 'ReadVariableNames', true);

% Assume the first three columns are the target log returns
dataArray = table2array(dataTable(:, 1:3));
variableNames = dataTable.Properties.VariableNames(1:3);

% Define horizons
horizons = [1, 3, 12]; % Monthly, Quarterly, Annual

% Initialize storage for results
summaryStats = struct();
monthlyCorrelationMatrix = [];
annualCorrelationMatrix = [];

% Compute summary statistics for each horizon
for h = horizons
    fprintf('Processing horizon: %d\n', h);
    
    % For log returns, no need to compute cumulative returns for horizon h
    logReturns = dataArray;
    
    % Prepare storage for statistics
    statsMean = nan(1, size(logReturns, 2));
    statsStd = nan(1, size(logReturns, 2));
    statsSkew = nan(1, size(logReturns, 2));
    statsKurt = nan(1, size(logReturns, 2));
    statsAR1 = nan(1, size(logReturns, 2));
    
    % Calculate statistics for each asset
    for i = 1:size(logReturns, 2)
        if h == 1
            assetReturns = logReturns(:, i);  % Monthly
        else
            assetReturns = computeCumulativeReturns(logReturns(:, i), h); % Quarterly, Annual
        end
        
        % Compute summary statistics
        statsMean(i) = mean(assetReturns);
        statsStd(i) = std(assetReturns);
        statsSkew(i) = skewness(assetReturns);
        statsKurt(i) = kurtosis(assetReturns);
        
        % Compute AR(1) coefficient using autocorrelation
        acf = autocorr(assetReturns, 1);
        statsAR1(i) = acf(2); % AR(1) coefficient
    end
    
    % Create a table for this horizon
    horizonTable = table(statsMean', statsStd', statsSkew', statsKurt', statsAR1', ...
        'VariableNames', {'Mean', 'Std', 'Skew', 'Kurt', 'AR1'}, ...
        'RowNames', strcat(variableNames, '_horizon_', num2str(h)));
    
    % Store results
    summaryStats.(['horizon_' num2str(h)]) = horizonTable;
end

% Compute the Monthly and Annual Correlation Matrices
monthlyCorrelationMatrix = corr(dataArray);
annualReturns = computeCumulativeReturns(dataArray, 12); % Annual returns
annualCorrelationMatrix = corr(annualReturns);

% Combine the matrices into a single matrix
numAssets = size(dataArray, 2);
combinedCorrelationMatrix = NaN(numAssets, numAssets);

% Fill the matrix with monthly correlations above the diagonal
combinedCorrelationMatrix(tril(true(numAssets), -1)) = NaN; % Lower triangle NaN
combinedCorrelationMatrix(triu(true(numAssets), 1)) = monthlyCorrelationMatrix(triu(true(numAssets), 1));

% Fill the matrix with annual correlations below the diagonal
combinedCorrelationMatrix(tril(true(numAssets), 0)) = annualCorrelationMatrix(tril(true(numAssets), 0));

% Create a table for the combined correlation matrix
combinedCorrelationTable = array2table(combinedCorrelationMatrix, ...
    'RowNames', variableNames, ...
    'VariableNames', variableNames);

% Display results
disp('Summary Statistics for Log Returns:');
disp('Panel A: Monthly');
disp(summaryStats.horizon_1);

disp('Panel B: Quarterly');
disp(summaryStats.horizon_3);

disp('Panel C: Annual');
disp(summaryStats.horizon_12);

disp('Panel D: Correlation Matrix');
disp(combinedCorrelationTable);

% Export Monthly Data and Plots
exportMonthlyDataAndPlots(dataArray, variableNames);

% Helper functions
function cumulativeReturns = computeCumulativeReturns(logReturns, h)
    % Compute cumulative returns for horizon h
    [T, numAssets] = size(logReturns);
    if T <= h
        error('Not enough data to compute returns for horizon %d', h);
    end
    
    cumulativeReturns = NaN(T - h, numAssets);
    for i = 1:numAssets
        for t = 1:(T - h)
            cumulativeReturns(t, i) = sum(logReturns(t:(t+h-1), i));
        end
    end
end

function exportMonthlyDataAndPlots(dataArray, variableNames)
    % Define indices for the three indices of interest
    indices = variableNames;
    numIndices = length(indices);
    
    % Create figures for log returns and volatility
    for i = 1:numIndices
        % Log returns are already given
        logReturns = dataArray(:, i);
        
        % Calculate rolling volatility with a window size of 12 months
        windowSize = 12;
        rollingVolatility = NaN(length(logReturns) - windowSize + 1, 1);
        for t = windowSize:length(logReturns)
            rollingVolatility(t - windowSize + 1) = std(logReturns(t - windowSize + 1:t));
        end
        
        % Plot log returns and rolling volatility
        figure;
        
        % Plot log returns
        subplot(2, 1, 1);
        plot(logReturns);
        title(['Monthly Log Returns of ', indices{i}]);
        xlabel('Time');
        ylabel('Log Return');
        grid on;
        
        % Plot rolling volatility
        subplot(2, 1, 2);
        timeVector = windowSize:length(logReturns);
        plot(timeVector, rollingVolatility);
        title(['Monthly Rolling Volatility of ', indices{i}]);
        xlabel('Time');
        ylabel('Volatility');
        grid on;
        
        % Save figures
        saveas(gcf, [indices{i}, '_Monthly_Plots.png']);
    end
    
    % Export the monthly data to a .csv file
    monthlyDataTable = array2table(dataArray, 'VariableNames', indices);
    writetable(monthlyDataTable, 'MonthlyData.csv');
end
