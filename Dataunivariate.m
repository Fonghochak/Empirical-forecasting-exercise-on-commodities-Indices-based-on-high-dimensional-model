% Define file path
filledDataPath = 'E:/OneDrive - Lancaster University/Msc MBF/2. ECO420 Dissertation/Week 3 Data/Final data/Data for matlab/dataMatrix_filled.xlsx';

% Check if the file exists
if ~isfile(filledDataPath)
    error('File does not exist. Please check the file path.');
end

% Read filled data
dataTable = readtable(filledDataPath, 'ReadVariableNames', true);

% Extract target and predictor data
targetDataArray = table2array(dataTable(:, 1:3)); % First three columns as target variables
predictorDataArray = table2array(dataTable(:, 4:end)); % Remaining columns as predictor variables

% Define horizons
horizons = [1, 3, 12]; % Monthly, Quarterly, Annual

% Define predictor groups (assuming predictors are in the same order as provided)
groups = {
    'Financial', {'logDP','logDY','logEP','logDE','BM','TBL','DFY','LTY','NTIS','INFL','LTR','TMS','DFR','SVAR'};
    'Group_1_Output_and_Income_16', {'RPI', 'W875RX1', 'INDPRO', 'IPFPNSS', 'IPFINAL', 'IPCONGD', 'IPDCONGD', 'IPNCONGD', 'IPBUSEQ', 'IPMAT', 'IPDMAT', 'IPNMAT', 'IPMANSICS', 'IPB51222S', 'IPFUELS', 'CUMFNS'};
    'Group_2_Labor_Market_31', {'HWI', 'HWIURATIO', 'CLF16OV', 'CE16OV', 'UNRATE', 'UEMPMEAN', 'UEMPLT5', 'UEMP5TO14', 'UEMP15OV', 'UEMP15T26', 'UEMP27OV', 'CLAIMSx', 'PAYEMS', 'USGOOD', 'CES1021000001', 'USCONS', 'MANEMP', 'DMANEMP', 'NDMANEMP', 'SRVPRD', 'USTPU', 'USWTRADE', 'USTRADE', 'USFIRE', 'USGOVT', 'CES0600000007', 'AWOTMAN', 'AWHMAN', 'CES0600000008', 'CES2000000008', 'CES3000000008'};
    'Group_3_Consumption_and_Orders_7', {'HOUSTS', 'HOUSTW', 'PERMIT', 'PERMITNE', 'PERMITMW', 'PERMITS', 'PERMITW'};
    'Group_4_Orders_and_Inventories_10', {'DPCERA3M086SBEA', 'CMRMTSPLx', 'RETAILx', 'ACOGNO', 'AMDMNOx', 'ANDENOx', 'AMDMUOx', 'BUSINVx', 'ISRATIOx', 'UMCSENTx'};
    'Group_5_Money_and_Credit_13', {'M1SL', 'M2SL', 'M2REAL', 'AMBSL', 'TOTRESNS', 'NONBORRES', 'BUSLOANS', 'REALLN', 'NONREVSL', 'CONSPI', 'DTCOLNVHFNM', 'DTCTHFNM', 'INVEST'};
    'Group_6_Interest_Rate_and_Exchange_Rates_22', {'FEDFUNDS', 'CP3Mx', 'TB3MS', 'TB6MS', 'GS1', 'GS5', 'GS10', 'AAA', 'BAA', 'COMPAPFFx', 'TB3SMFFM', 'TB6SMFFM', 'T1YFFM', 'T5YFFM', 'T10YFFM', 'AAAFFM', 'BAAFFM', 'TWEXAFEGSMTHx', 'EXSZUSx', 'EXJPUSx', 'EXUSUKx', 'EXCAUSx'};
    'Group_7_Prices_20', {'WPSFD49207', 'WPSFD49502', 'WPSID61', 'WPSID62', 'OILPRICEx', 'PPICMM', 'CPIAUCSL', 'CPIAPPSL', 'CPITRNSL', 'CPIMEDSL', 'CUSR0000SAC', 'CUSR0000SAD', 'CUSR0000SAS', 'CPIULFSL', 'CUSR0000SA0L2', 'CUSR0000SA0L5', 'PCEPI', 'DDURRG3M086SBEA', 'DNDGRG3M086SBEA', 'DSERRG3M086SBEA'};
    'Group_8_Stock_Market_3', {'SP500', 'SPDY', 'SPPE'};
};

% Initialize storage for results
results = [];

% Compute summary statistics and regression results for each horizon
for h = horizons
    fprintf('Processing horizon: %d\n', h);
    
    % Compute cumulative returns for horizon h
    cumulativeReturns = computeCumulativeReturns(targetDataArray, h);
    
    % Ensure predictor data matches the number of rows in cumulative returns
    numObs = min(size(cumulativeReturns, 1), size(predictorDataArray, 1));
    cumulativeReturns = cumulativeReturns(1:numObs, :);
    predictorDataArray = predictorDataArray(1:numObs, :);
    
    % Run regressions and compute significance for each group of predictors
    for g = 1:size(groups, 1)
        groupName = groups{g, 1};
        predictors = groups{g, 2};
        
        for p = 1:length(predictors)
            predictorName = predictors{p};
            predictorIndex = find(strcmp(dataTable.Properties.VariableNames, predictorName));
            if isempty(predictorIndex)
                error('Predictor %s not found in the data table.', predictorName);
            end
            predictorData = predictorDataArray(:, predictorIndex - 3); % Adjust index for the predictorDataArray
            
            % Run OLS regression for each asset
            for assetIdx = 1:3
                assetName = dataTable.Properties.VariableNames{assetIdx};
                assetReturns = cumulativeReturns(:, assetIdx);
                
                % Ensure predictor data and asset returns match in length
                if length(assetReturns) ~= length(predictorData)
                    error('Mismatch in number of observations between asset returns and predictor data.');
                end
                
                % Run OLS regression
                [b, ~, ~, ~, stats] = regress(assetReturns, [ones(size(predictorData, 1), 1), predictorData]);
                
                % Compute bootstrap p-values
                pVal = bootstrapPValue(assetReturns, predictorData, b, 1000);
                
                % Store results
                results = [results; {h, groupName, predictorName, assetName, b(2), pVal}];
            end
        end
    end
end

% Convert results to table and save to CSV
resultsTable = cell2table(results, 'VariableNames', {'Horizon', 'Group', 'Predictor', 'Target', 'Beta', 'pValue'});
writetable(resultsTable, 'univariate_volatility_results.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions

function cumulativeReturns = computeCumulativeReturns(prices, h)
    % Compute cumulative returns for horizon h
    [T, numAssets] = size(prices);
    if T <= h
        error('Not enough data to compute returns for horizon %d', h);
    end
    
    cumulativeReturns = NaN(T - h, numAssets);
    for i = 1:numAssets
        for t = 1:(T - h)
            cumulativeReturns(t, i) = (prices(t + h, i) - prices(t, i)) / prices(t, i);
        end
    end
end

function pVal = bootstrapPValue(cumulativeReturns, predictorData, beta, numBootstrap)
    % Compute p-values using bootstrap method
    numSamples = size(cumulativeReturns, 1);
    bootstrapBetas = zeros(numBootstrap, 1);
    
    for b = 1:numBootstrap
        % Resample with replacement
        indices = randsample(numSamples, numSamples, true);
        resampledReturns = cumulativeReturns(indices);
        resampledPredictor = predictorData(indices);
        
        % Run OLS regression on resampled data
        [bResampled, ~, ~, ~, ~] = regress(resampledReturns, [ones(size(resampledPredictor, 1), 1), resampledPredictor]);
        bootstrapBetas(b) = bResampled(2);
    end
    
    % Calculate p-value based on bootstrap distribution
    pVal = mean(abs(bootstrapBetas) >= abs(beta(2)));
end
