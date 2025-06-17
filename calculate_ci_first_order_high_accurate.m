function [CI, UB] = calculate_ci_first_order_high_accurate(R_hat, params, n, confidence_level)
% =========================================================================
% 使用一阶逼近法计算CEP的置信区间和置信上界的高精度计算。慢！
% =========================================================================
%
%   - 严格遵循 GJB 6289-2008 附录B 中的解析公式计算偏导数，
%     摒弃数值差分方法，以获得最高精度。
%
%--------------------------------------------------------------------------

    mu1 = params.mu1;
    mu2 = params.mu2;
    s1 = params.s1;
    s2 = params.s2;
    s1_sq = s1^2;
    s2_sq = s2^2;
    
    % --- 1. 计算附录B中的常数和中间变量 ---
    c_exp_term = exp(-0.5 * (mu1^2/s1_sq + mu2^2/s2_sq));
    c = (1/(s1*s2)) * c_exp_term;
    
    % c_i (公式 B.6, B.7)
    c1 = -(mu1 / s1_sq) * c;
    c2 = -(mu2 / s2_sq) * c;
    c3 = (c / (2*s1_sq)) * ((mu1^2/s1_sq) - 1);
    c4 = (c / (2*s2_sq)) * ((mu2^2/s2_sq) - 1);

    % T 和 T(i) (公式 B.5, B.8 - B.11)
    % 需要计算相关的 L 积分项
    fprintf('    正在计算附录B中的积分项...\n');
    L1_00 = calculate_L_integrals(1, [0, 0], params, R_hat);
    L2_210 = calculate_L_integrals(2, [2, 1, 0], params, R_hat);
    L2_201 = calculate_L_integrals(2, [2, 0, 1], params, R_hat);
    L2_320 = calculate_L_integrals(2, [3, 2, 0], params, R_hat);
    L2_302 = calculate_L_integrals(2, [3, 0, 2], params, R_hat);
    
    T = L1_00 * c * R_hat;
    T1 = L2_210 / s1_sq;
    T2 = L2_201 / s2_sq;
    T3 = (1 / (2*s1_sq^2)) * (L2_320 - 2*mu1*L2_210);
    T4 = (1 / (2*s2_sq^2)) * (L2_302 - 2*mu2*L2_201);
    fprintf('    积分项计算完成。\n');

    % --- 2. 计算一阶偏导数 R_i = dR/d(theta_i) (公式 B.1) ---
    % 注意：规范中的 theta_i 是 [mu1, mu2, s1^2, s2^2]
    
    % dR/d(mu1)
    dR_dmu1 = -(c1/(2*c) + c * T1) / T;
    
    % dR/d(mu2)
    dR_dmu2 = -(c2/(2*c) + c * T2) / T;
    
    % dR/d(s1^2)
    dR_ds1sq = -(c3/(2*c) + c * T3) / T;
    
    % dR/d(s2^2)
    dR_ds2sq = -(c4/(2*c) + c * T4) / T;

    % --- 3. 计算 H_psi ---
    % 变换函数 psi(x) = ln(x), 其导数为 1/x
    psi_prime_R = 1 / R_hat;

    % 计算 psi_i
    psi_1 = psi_prime_R * dR_dmu1;
    psi_2 = psi_prime_R * dR_dmu2;
    psi_3 = psi_prime_R * dR_ds1sq;
    psi_4 = psi_prime_R * dR_ds2sq;
    
    % 规范中的 v_i
    v1 = s1_sq;
    v2 = s2_sq;
    v3 = 2 * (s1_sq^2);
    v4 = 2 * (s2_sq^2);
    
    % 公式(24)
    H_psi_sq_sum = v1 * psi_1^2 + v2 * psi_2^2 + v3 * psi_3^2 + v4 * psi_4^2;
    H_psi = 1./sqrt(H_psi_sq_sum);
    
    % --- 4. 计算置信区间 ---
    alpha = 1 - confidence_level;
    z_score_lower = norminv(alpha / 2);
    z_score_upper = norminv(1 - alpha / 2);
    
    factor_lower = z_score_lower / sqrt(n * H_psi^2);
    factor_upper = z_score_upper / sqrt(n * H_psi^2);
    
    % 公式(28, 29)
    CI_lower = R_hat * exp(factor_lower);
    CI_upper = R_hat * exp(factor_upper);
    CI = [CI_lower, CI_upper];
    
    % 公式(30)
    z_score_ub = norminv(1 - alpha);
    factor_ub = z_score_ub / sqrt(n * H_psi^2);
    UB = R_hat * exp(factor_ub);
end