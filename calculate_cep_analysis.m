function results = calculate_cep_analysis(x, z, confidence_level, M_boot, test_alpha)
% =========================================================================
% GJB 6289-2008 核心分析函数
% =========================================================================
%
% 功能:
%   对输入数据进行正态性检验 (Shapiro-Wilk) 和独立性检验 (Pearson)。
%   根据输入的落点偏差数据(x,z)，完成 GJB 6289-2008 中定义的
%   CEP点估计和区间估计。
% 输入:
%   x                - 纵向落点偏差 (行向量)
%   z                - 横向落点偏差 (行向量)
%   confidence_level - 置信水平 (e.g., 0.95)
%   M_boot           - Bootstrap 模拟次数 (e.g., 2000)
%   test_alpha       - 假设检验的显著性水平 (e.g., 0.10)
%
% 输出:
%   results          - 包含所有计算结果的结构体
%
%--------------------------------------------------------------------------

    % --- 0. 数据检验 (GJB 4.3节) ---
    fprintf('  正在进行数据有效性检验 (显著性水平 α = %.2f)...\n', test_alpha);
    
    % a) 正态性检验 (Shapiro-Wilk Test)
    [h_norm_x, p_norm_x] = sf_test(x, test_alpha);
    [h_norm_z, p_norm_z] = sf_test(z, test_alpha);
    
    results.tests.normality.pass_x = ~h_norm_x;
    results.tests.normality.p_value_x = p_norm_x;
    results.tests.normality.pass_z = ~h_norm_z;
    results.tests.normality.p_value_z = p_norm_z;
    
    if h_norm_x || h_norm_z
        fprintf('  [检验失败] 数据未通过正态性检验。\n');
        if h_norm_x, fprintf('    - 纵向(X)数据 P-value = %.4f < %.2f\n', p_norm_x, test_alpha); end
        if h_norm_z, fprintf('    - 横向(Z)数据 P-value = %.4f < %.2f\n', p_norm_z, test_alpha); end
        fprintf('  根据GJB 6289-2008，后续计算不适用。请检查数据源。\n');
        results.error = '数据未通过正态性检验';
        return;
    else
        fprintf('  [检验通过] 数据通过正态性检验。\n');
    end

    % b) 独立性检验 (Pearson Correlation)
    [rho_matrix, p_corr_matrix] = corrcoef(x, z);
    rho = rho_matrix(1, 2);
    p_corr = p_corr_matrix(1, 2);
    
    is_independent = (p_corr >= test_alpha);
    results.tests.independence.is_independent = is_independent;
    results.tests.independence.p_value = p_corr;
    
    if is_independent
        fprintf('  [检验通过] 数据通过独立性检验 (P-value = %.4f >= %.2f)。\n\n', p_corr, test_alpha);
        rho = 0; % 强制为0，因为统计上不显著相关
    else
        fprintf('  [检验表明不独立] 数据未通过独立性检验 (P-value = %.4f < %.2f)。\n', p_corr, test_alpha);
        fprintf('  将采用相关模型进行正交变换。\n\n');
    end

    % --- 1. 计算基础统计参数 (对应规范第6章) ---
    n = length(x);
    mux = mean(x);
    muz = mean(z);
    sx = std(x);
    sz = std(z);
    
    results.basic_stats.n = n;
    results.basic_stats.mux = mux;
    results.basic_stats.muz = muz;
    results.basic_stats.sx = sx;
    results.basic_stats.sz = sz;
    results.basic_stats.rho = rho;

    % --- 2. 正交变换 (若需要) ---
    % 根据规范第5.3节, 当相关性检验不通过时时，进行变换
    if ~is_independent
        % 公式(14)
        xi = 0.5 * atan((2 * rho * sx * sz) / (sx^2 - sz^2));
        
        % 公式(15)
        mu_x_prime = mux * cos(xi) + muz * sin(xi);
        mu_z_prime = muz * cos(xi) - mux * sin(xi);
        
        factor = sqrt(1 - rho^2);
        term_s1_sq = sz^2 * cos(xi)^2 + sx^2 * sin(xi)^2 - rho * sx * sz * sin(2*xi);
        term_s2_sq = sz^2 * sin(xi)^2 + sx^2 * cos(xi)^2 + rho * sx * sz * sin(2*xi);
        
        s_x_prime = sx * sz * factor / sqrt(term_s1_sq);
        s_z_prime = sx * sz * factor / sqrt(term_s2_sq);

        % 统一参数 (公式 16-19)
        mu1 = mu_x_prime;
        mu2 = mu_z_prime;
        s1 = s_x_prime;
        s2 = s_z_prime;
    else
        % 独立情况
        mu1 = mux;
        mu2 = muz;
        s1 = sx;
        s2 = sz;
    end
    
    % 将用于计算的参数存入结构体
    params.mu1 = mu1;
    params.mu2 = mu2;
    params.s1 = s1;
    params.s2 = s2;
    results.params = params;

    % --- 3. CEP 点估计 ---
    R_hat = calculate_cep_plugin(params);
    results.point_estimate.R_hat = R_hat;
    
    % --- 4. CEP 区间估计 ---
    alpha = 1 - confidence_level;
    
    % 方法1: 参数自助法
    fprintf('  正在进行参数自助法计算 (模拟 %d 次)...\n', M_boot);
    [CI_boot, UB_boot] = calculate_ci_bootstrap(R_hat, params, n, confidence_level, M_boot);
    results.ci_bootstrap.CI = CI_boot;
    results.ci_bootstrap.UB = UB_boot;
    fprintf('  参数自助法计算完成。\n');

    % 方法2: 一阶逼近法
    fprintf('  正在进行一阶逼近法计算...\n');
    [CI_approx, UB_approx] = calculate_ci_first_order(R_hat, params, n, confidence_level);
    results.ci_first_order.CI = CI_approx;
    results.ci_first_order.UB = UB_approx;
    fprintf('  一阶逼近法计算完成。\n\n');
end