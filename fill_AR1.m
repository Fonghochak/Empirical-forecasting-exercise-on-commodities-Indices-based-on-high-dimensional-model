function filled_series = fill_AR1(series, phi, sigma)
    % 填补时间序列中的缺失值使用 AR(1) 模型
    % 输入:
    %   series - 包含缺失值的时间序列
    %   phi - AR(1) 系列的系数
    %   sigma - 噪声的标准差
    % 输出:
    %   filled_series - 填补缺失值后的时间序列

    n = length(series);
    for t = 2:n
        if isnan(series(t))
            series(t) = phi * series(t-1) + normrnd(0, sigma);
        end
    end
    filled_series = series;
end
