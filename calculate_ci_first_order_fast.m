function [CI, UB] = calculate_ci_first_order_fast(R_hat, params, n, confidence_level)
% =========================================================================
% 使用一阶逼近法计算CEP的置信区间和置信上界
% =========================================================================
%
% 涉及公式: GJB 6289-2008, 第7.3节, 公式(28, 29, 30)
%
%--------------------------------------------------------------------------
    mu1 = params.mu1;
    mu2 = params.mu2;
    s1 = params.s1;
    s2 = params.s2;
    
    % --- 1. 计算 H_psi ---
    % 规范中的 theta = [mu1, mu2, s1^2, s2^2]
    
    % 使用数值微分计算 R 对 theta 的偏导数，这里没有直接使用附录B的复杂公式
    h = 1e-6; % 微小扰动

    % dR/d(mu1)
    params_h = params; params_h.mu1 = mu1 + h;
    R_h_mu1 = calculate_cep_plugin(params_h);
    dR_dmu1 = (R_h_mu1 - R_hat) / h;

    % dR/d(mu2)
    params_h = params; params_h.mu2 = mu2 + h;
    R_h_mu2 = calculate_cep_plugin(params_h);
    dR_dmu2 = (R_h_mu2 - R_hat) / h;
    
    % dR/d(s1^2)
    params_h = params; params_h.s1 = sqrt(s1^2 + h);
    R_h_s1sq = calculate_cep_plugin(params_h);
    dR_ds1sq = (R_h_s1sq - R_hat) / h;

    % dR/d(s2^2)
    params_h = params; params_h.s2 = sqrt(s2^2 + h);
    R_h_s2sq = calculate_cep_plugin(params_h);
    dR_ds2sq = (R_h_s2sq - R_hat) / h;

    % 变换函数 psi(x) = ln(x), 其导数为 1/x
    psi_prime_R = 1 / R_hat;

    % 计算 psi_i
    psi_1 = psi_prime_R * dR_dmu1;
    psi_2 = psi_prime_R * dR_dmu2;
    psi_3 = psi_prime_R * dR_ds1sq;
    psi_4 = psi_prime_R * dR_ds2sq;
    
    % 规范中的 v_i
    v1 = s1^2;
    v2 = s2^2;
    v3 = 2 * (s1^4);
    v4 = 2 * (s2^4);
    
    % 公式(24)
    H_psi_sq = v1 * psi_1^2 + v2 * psi_2^2 + v3 * psi_3^2 + v4 * psi_4^2;
    H_psi = 1./sqrt(H_psi_sq);

    % --- 2. 计算置信区间 ---
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