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

% Define target names for better visualization
target_names = {'SPGSEN Index', 'SPGSIN Index', 'SPGSAGS Index'};

% Initialize result storage for AR(1) model evaluation metrics
num_targets = size(y_targets, 2);
rmse_values_ar1 = zeros(num_targets, 1);
mae_values_ar1 = zeros(num_targets, 1);
mape_values_ar1 = zeros(num_targets, 1);

for i = 1:num_targets
    y_it_in_sample = y_in_sample(:, i);
    y_it_out_sample = y_out_sample(:, i);
    
    % Fit AR(1) model to in-sample data
    ar1_mdl = arima(1, 0, 0);
    ar1_fit = estimate(ar1_mdl, y_it_in_sample);
    
    % Use the AR(1) model to predict out-sample data
    ar1_predictions = forecast(ar1_fit, out_sample_obs, 'Y0', y_it_in_sample);
    
    % Calculate residuals for AR(1) model
    residuals_ar1 = y_it_out_sample - ar1_predictions;
    
    % Calculate AR(1) model evaluation metrics
    rmse_values_ar1(i) = sqrt(mean(residuals_ar1.^2));
    mae_values_ar1(i) = mean(abs(residuals_ar1));
    mape_values_ar1(i) = mean(abs(residuals_ar1 ./ y_it_out_sample) * 100);
    
    % Display AR(1) results
    fprintf('AR(1) Model for Target %d (%s):\n', i, target_names{i});
    fprintf('  AR(1) RMSE = %.4f\n', rmse_values_ar1(i));
    fprintf('  AR(1) MAE = %.4f\n', mae_values_ar1(i));
    fprintf('  AR(1) MAPE = %.4f\n', mape_values_ar1(i));
    
    % Plot original vs predicted data
    figure;
    hold on;
    plot(y_it_out_sample, 'b', 'DisplayName', 'Original Data');
    plot(ar1_predictions, 'g--', 'DisplayName', 'AR(1) Predicted Data');
    hold off;
    
    % Set plot title and labels
    title(['Original vs AR(1) Predicted Data for ', target_names{i}]);
    xlabel('Time');
    ylabel('Value');
    legend('show');
    
    % Save the figure
    saveas(gcf, ['AR1_Prediction_vs_Original_', target_names{i}, '.png']);
end
