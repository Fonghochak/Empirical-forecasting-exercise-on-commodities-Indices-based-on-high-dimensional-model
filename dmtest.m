% Define the Diebold-Mariano test function
function dm_stat = dmtest(e1, e2, h)
    d = e1 - e2;
    mean_d = mean(d);
    var_d = var(d);
    if var_d == 0
        dm_stat = NaN;
    else
        dm_stat = mean_d / sqrt(var_d / length(d));
    end
end
