function [icp1, icp2, icp3] = calculate_ic(Xfactor, maxLag)
    [T, N] = size(Xfactor);
    icp1 = NaN(maxLag, N);
    icp2 = NaN(maxLag, N);
    icp3 = NaN(maxLag, N);

    for lag = 1:maxLag
        % 假设使用固定滞后量计算信息准则
        % 你可以使用模型来计算实际的误差项，这里仅演示格式
        [F, ~, ~] = tmp(Xfactor, lag);
        
        % 计算残差平方和作为误差
        residuals = Xfactor - F; % 示例，实际中可能需要使用具体模型计算残差
        err_sum = sum(residuals(:).^2);

        % 避免对 0 取对数
        if err_sum <= 0
            log_err_sum = -Inf;
        else
            log_err_sum = log(err_sum);
        end
        
        % 这里的罚项需要根据实际模型进行调整
        pf1 = ((N + T) / (N * T)) * log((N * T) / (N + T));
        pf2 = ((N + T) / (N * T)) * log(min([N, T]));
        pf3 = log(min([N, T])) / min([N, T]);

        % 计算信息准则
        icp1(lag, :) = log_err_sum + lag * pf1;
        icp2(lag, :) = log_err_sum + lag * pf2;
        icp3(lag, :) = log_err_sum + lag * pf3;
    end
end
