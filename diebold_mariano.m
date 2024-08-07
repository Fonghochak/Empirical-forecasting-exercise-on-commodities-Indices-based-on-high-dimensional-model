function pval = diebold_mariano(y_true, forecast1, forecast2)
    % Diebold-Mariano test for comparing predictive accuracy of two models
    % Inputs:
    % y_true: True values
    % forecast1: Forecast from benchmark model (historical average)
    % forecast2: Forecast from the complex model (PCA-based)

    % Calculate forecast errors
    e1 = y_true - forecast1;
    e2 = y_true - forecast2;
    
    % Calculate the loss differential
    d = e1.^2 - e2.^2;
    
    % Compute test statistic
    T = length(y_true);
    DM_stat = mean(d) / (std(d) / sqrt(T));
    
    % Compute p-value
    pval = 2 * (1 - normcdf(abs(DM_stat)));
end
