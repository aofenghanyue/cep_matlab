function [CI, UB] = calculate_ci_first_order(R_hat, params, n, confidence_level)
% =========================================================================
%   method = 1: 数值差分法 (速度快，精度较低)。
%   method = 2: 解析公式法 (GJB附录B，速度慢，精度高)。
method = 2;
if method == 1
    [CI, UB] = calculate_ci_first_order_fast(R_hat, params, n, confidence_level);
else
    [CI, UB] = calculate_ci_first_order_high_accurate(R_hat, params, n, confidence_level);
end
end