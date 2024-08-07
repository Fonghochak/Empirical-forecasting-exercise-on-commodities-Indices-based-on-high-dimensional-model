function pval = giacomini_white(y_true, forecast1, forecast2, X_out_sample)
    % Giacomini-White test for comparing predictive accuracy of two models
    % Inputs:
    % y_true: True values
    % forecast1: Forecast from benchmark model (historical average)
    % forecast2: Forecast from the complex model (PCA-based)
    % X_out_sample: Out-sample predictor variables
    
    % Calculate forecast errors
    e1 = y_true - forecast1;
    e2 = y_true - forecast2;
    
    % Calculate the loss differential
    d = e1.^2 - e2.^2;
    
    % Regress the loss differential on a constant and out-sample predictors
    T = length(y_true);
    Z = [ones(T, 1), X_out_sample];
    beta = (Z' * Z) \ (Z' * d);
    residuals = d - Z * beta;
    
    % Compute test statistic
    GW_stat = beta' * (Z' * Z / T) * beta / (residuals' * residuals / (T - size(Z, 2)));
    
    % Compute p-value
    pval = 1 - chi2cdf(GW_stat, size(Z, 2));
end
