function [F, coeff, explained] = tmp(X, K)
    % 确保输入 X 是数值矩阵，没有 NaN 或 Inf 值
    if any(isnan(X(:))) || any(isinf(X(:)))
        error('输入数据包含 NaN 或 Inf 值');
    end

    % 执行主成分分析
    [coeff, ~, latent] = pca(X, 'NumComponents', K);
    F = X * coeff; % 计算因子得分
    explained = 100 * latent / sum(latent); % 解释的方差比例
end
