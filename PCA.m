%%
%panety function
% Clear workspace
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

% Define the groups based on the department
groups = {
    'Financial', {'logDP','logDY','logEP','logDE','BM','TBL','DFY','LTY','NTIS','INFL','LTR','TMS','DFR','SVAR'};
    'Group 1: Output and Income 16', {'RPI', 'W875RX1', 'INDPRO', 'IPFPNSS', 'IPFINAL', 'IPCONGD', 'IPDCONGD', 'IPNCONGD', 'IPBUSEQ', 'IPMAT', 'IPDMAT', 'IPNMAT', 'IPMANSICS', 'IPB51222S', 'IPFUELS', 'CUMFNS'};
    'Group 2: Labor Market 31', {'HWI', 'HWIURATIO', 'CLF16OV', 'CE16OV', 'UNRATE', 'UEMPMEAN', 'UEMPLT5', 'UEMP5TO14', 'UEMP15OV', 'UEMP15T26', 'UEMP27OV', 'CLAIMSx', 'PAYEMS', 'USGOOD', 'CES1021000001', 'USCONS', 'MANEMP', 'DMANEMP', 'NDMANEMP', 'SRVPRD', 'USTPU', 'USWTRADE', 'USTRADE', 'USFIRE', 'USGOVT', 'CES0600000007', 'AWOTMAN', 'AWHMAN', 'CES0600000008', 'CES2000000008', 'CES3000000008'};
    'Group 3: Consumption and Orders 7', {'HOUSTS', 'HOUSTW', 'PERMIT', 'PERMITNE', 'PERMITMW', 'PERMITS', 'PERMITW'};
    'Group 4: Orders and Inventories 10', {'DPCERA3M086SBEA', 'CMRMTSPLx', 'RETAILx', 'ACOGNO', 'AMDMNOx', 'ANDENOx', 'AMDMUOx', 'BUSINVx', 'ISRATIOx', 'UMCSENTx'};
    'Group 5: Money and Credit 13', {'M1SL', 'M2SL', 'M2REAL', 'AMBSL', 'TOTRESNS', 'NONBORRES', 'BUSLOANS', 'REALLN', 'NONREVSL', 'CONSPI', 'DTCOLNVHFNM', 'DTCTHFNM', 'INVEST'};
    'Group 6: Interest Rate and Exchange Rates 22', {'FEDFUNDS', 'CP3Mx', 'TB3MS', 'TB6MS', 'GS1', 'GS5', 'GS10', 'AAA', 'BAA', 'COMPAPFFx', 'TB3SMFFM', 'TB6SMFFM', 'T1YFFM', 'T5YFFM', 'T10YFFM', 'AAAFFM', 'BAAFFM', 'TWEXAFEGSMTHx', 'EXSZUSx', 'EXJPUSx', 'EXUSUKx', 'EXCAUSx'};
    'Group 7: Prices 20', {'WPSFD49207', 'WPSFD49502', 'WPSID61', 'WPSID62', 'OILPRICEx', 'PPICMM', 'CPIAUCSL', 'CPIAPPSL', 'CPITRNSL', 'CPIMEDSL', 'CUSR0000SAC', 'CUSR0000SAD', 'CUSR0000SAS', 'CPIULFSL', 'CUSR0000SA0L2', 'CUSR0000SA0L5', 'PCEPI', 'DDURRG3M086SBEA', 'DNDGRG3M086SBEA', 'DSERRG3M086SBEA'};
    'Group 8: Stock Market 3', {'SP500', 'SPDY', 'SPPE'};
    %'VIX Index', {'VIXCLSx'}
};

% Initialize result storage
all_r_squared_tables = cell(size(groups, 1), 1);

% Define the penalty function g(N, T)
[N, T] = size(dataMatrix);
g = @(N, T) (N + T) * min(N, T) / (N * T);

% Loop over each group
for g_idx = 1:size(groups, 1)
    group_name = groups{g_idx, 1};
    group_vars = groups{g_idx, 2};
    
    % Find the indices of the current group variables in the data table
    [found, group_indices] = ismember(group_vars, dataTable.Properties.VariableNames);
    
    % Check for missing variables
    if any(~found)
        missing_vars = group_vars(~found);
        warning('The following variables are missing in the data table for %s:\n%s', group_name, strjoin(missing_vars, ', '));
        % Remove missing variables from the group
        group_vars = group_vars(found);
        group_indices = group_indices(found);
    end
    
    % Extract predictor variables for the current group
    X_group = dataMatrix(:, group_indices);
    
    % Replace infinite values with NaN
    X_group(isinf(X_group)) = NaN;
    
    % Ensure X_group is not empty and has more than one column
    if isempty(X_group) || size(X_group, 2) < 2 || all(isnan(X_group), 'all')
        warning('Input matrix X_group is empty, has less than two columns, or contains only NaNs for %s. Skipping this group.', group_name);
        continue; % Skip this group and move to the next one
    end
    
    % Perform PCA to find latent factors
    try
        [coeff, score, latent, tsquared, explained] = pca(X_group, 'Rows', 'pairwise'); % PCA as a factor extraction method
    catch ME
        warning('PCA failed for %s: %s', group_name, ME.message);
        continue; % Skip this group and move to the next one
    end
    
    % Calculate the information criterion for each possible number of factors
    IC = zeros(size(score, 2), 1);
    max_factors = min(size(X_group, 2), size(X_group, 1) - 1); % Ensure the maximum number of factors does not exceed the dimensions of X_group
    for r = 1:max_factors
        F_hat = score(:, 1:r);
        Lambda_hat = coeff(:, 1:r);
        V_r = sum(sum((X_group - F_hat * Lambda_hat').^2));
        IC(r) = log(V_r) + r * g(N, T);
    end
    
    % Find the optimal number of factors
    [~, optimal_r] = min(IC(1:max_factors));
    
    % Use the optimal number of factors as the latent factors
    latent_factor = score(:, 1:optimal_r);
    
    % Initialize R-squared values storage
    r_squared_values = nan(1, size(y_targets, 2)); % One latent factor for each target variable
    
    % Perform regression analysis
    for i = 1:size(y_targets, 2)
        y_it = y_targets(:, i);
        
        % Perform ordinary least squares regression
        mdl = fitlm(latent_factor, y_it);
        
        % Calculate R-squared value
        r_squared_values(i) = mdl.Rsquared.Ordinary;
    end
    
    % Create results table for the current group
    r_squared_table = array2table(r_squared_values, 'VariableNames', {'SPGSEN Index', 'SPGSIN Index', 'SPGSAGS Index'}, 'RowNames', {group_name});
    
    % Store the table in the results cell array
    all_r_squared_tables{g_idx} = r_squared_table;
    
    % Display the results
    disp(['Table for ', group_name]);
    disp(r_squared_table);
    
end

% Combine all results into a single table for easier review
combined_r_squared_table = vertcat(all_r_squared_tables{:});

% Display combined results
disp('Combined R-squared Table');
disp(combined_r_squared_table);

% Save combined results to CSV
writetable(combined_r_squared_table, 'Combined_R_squared_table.csv', 'WriteRowNames', true);

% Prepare data for heatmap
r_squared_matrix = table2array(combined_r_squared_table);
row_names = combined_r_squared_table.Properties.RowNames;
col_names = combined_r_squared_table.Properties.VariableNames;

% Create heatmap
figure;
h = heatmap(r_squared_matrix, 'Colormap', parula, 'ColorbarVisible', 'on'); % Use 'parula' for a gradient colormap
h.XDisplayLabels = col_names;
h.YDisplayLabels = row_names;
h.Title = 'R-squared Heatmap';
h.XLabel = 'Targets';
h.YLabel = 'Groups';

% Save heatmap as image
saveas(gcf, 'R_squared_heatmap.png');
%%
% 清理工作空间和命令窗口
clear;
clc;

% 定义文件路径
filledDataPath = 'E:/OneDrive - Lancaster University/Msc MBF/2. ECO420 Dissertation/Week 3 Data/Final data/Data for matlab/dataMatrix_filled.xlsx';

% 检查文件是否存在
if ~isfile(filledDataPath)
    error('File does not exist. Please check the file path.');
end

% 读取数据，确保原始列名被保留
dataTable = readtable(filledDataPath, 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');
dataMatrix = readmatrix(filledDataPath);

% 数据标准化
dataMatrix = zscore(dataMatrix);

% 假设前三列是目标变量
y_targets = dataMatrix(:, 1:3);

% 剩余列假设为预测变量
predictors = dataMatrix(:, 4:end);
[N, T] = size(predictors);

% 对包含无穷值的数据进行处理
predictors(isinf(predictors)) = NaN;

% 执行PCA
[coeff, score, latent] = pca(predictors, 'Algorithm', 'eig', 'Rows', 'complete');

% 初始化R^2存储矩阵
total_r_squared = zeros(size(score, 2), size(predictors, 2));

% 计算每个因子对每个预测变量的R^2
for i = 1:size(score, 2)
    for j = 1:size(predictors, 2)
        mdl = fitlm(score(:, i), predictors(:, j));
        total_r_squared(i, j) = mdl.Rsquared.Ordinary;
    end
end

% 选择R^2最大的前五个因子
[~, sorted_indices] = sort(sum(total_r_squared, 2), 'descend');
top_factors = sorted_indices(1:min(5, length(sorted_indices)));
latent_factors = score(:, top_factors);

% 计算并显示每个因子对应的变量和其R^2
for i = 1:length(top_factors)
    fprintf('Latent Factor %d:\n', i);
    for j = 1:size(predictors, 2)
        if total_r_squared(top_factors(i), j) > 0.5  % 根据需要调整阈值
            fprintf('\t%s: R^2 = %.2f\n', dataTable.Properties.VariableNames{3 + j}, total_r_squared(top_factors(i), j));
        end
    end
end

% 对每个目标序列使用潜在因子进行回归分析
r_squared_values = zeros(size(y_targets, 2), 1);
for i = 1:size(y_targets, 2)
    mdl = fitlm(latent_factors, y_targets(:, i));
    r_squared_values(i) = mdl.Rsquared.Ordinary;
end

% 输出回归分析的R-squared值
target_names = dataTable.Properties.VariableNames(1:3);
r_squared_table = array2table(r_squared_values, 'VariableNames', target_names, 'RowNames', {'Combined Factors'});
disp('R-squared values for each target variable using all factors:');
disp(r_squared_table);


