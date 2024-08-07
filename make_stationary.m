% Function to check stationarity and difference if needed
function [data, diff_order] = make_stationary(data)
    diff_order = 0;
    p_value_threshold = 0.05;
    while true
        % Perform ADF test
        [~, p_value, ~, ~] = adftest(data);
        if p_value < p_value_threshold
            break;
        else
            data = diff(data);
            diff_order = diff_order + 1;
        end
    end
end