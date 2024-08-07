function pval = clarkwest(y_true, forecast1, forecast2, h)
    % Clark-West test for comparing predictive accuracy of two models
    % Inputs:
    % y_true: True values
    % forecast1: Forecast from benchmark model (historical average)
    % forecast2: Forecast from the complex model (PCA-based)
    % h: Forecast horizon (set to 1 for one-step ahead forecast)
    
    % Calculate forecast errors
    e1 = y_true - forecast1;
    e2 = y_true - forecast2;
    
    % Calculate adjusted forecast error for complex model
    f_tilde = e2.^2 - (e1.^2 - e2.^2);
    
    % Compute test statistic
    T = length(y_true);
    CW_stat = mean(f_tilde) / (std(f_tilde) / sqrt(T));
    
    % Compute p-value
    pval = 1 - normcdf(CW_stat);
end
