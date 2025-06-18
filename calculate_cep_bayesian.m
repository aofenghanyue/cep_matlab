function results = calculate_cep_bayesian(x_test, z_test, x_supp, z_supp, confidence_level, M)
% =========================================================================
% 使用贝叶斯参数自助法计算CEP的点估计和置信区间 (GJB 6289-2008, 第8章)
% =========================================================================
%
% 功能:
%   适用于小样本情况，结合补充样本（先验信息）和试验样本，
%   给出更稳健的CEP评定结果。
%
% 输入:
%   x_test, z_test   - 试验样本数据 (行向量, n <= 3)
%   x_supp, z_supp   - 补充样本数据 (行向量, 即先验信息)
%   confidence_level - 置信水平 (e.g., 0.80)
%   M                - 模拟次数 (e.g., 2000)
%
% 输出:
%   results          - 包含贝叶斯点估计、置信区间和上界的结构体
%
%--------------------------------------------------------------------------

    % --- 1. 计算后验分布的参数 ---
    % 我们将为 x 和 z 方向分别计算后验参数
    post_params_x = get_posterior_params(x_test, x_supp);
    post_params_z = get_posterior_params(z_test, z_supp);

    % --- 2. 从后验分布中进行蒙特卡洛抽样 ---
    fprintf('  正在进行贝叶斯参数抽样 (模拟 %d 次)...\n', M);
    
    R_bpt_samples = zeros(M, 1);
    
    % 创建并行池
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool;
    end
    
    % 为了在 parfor 中使用，将结构体中的字段提取为独立变量
    % 这些将被广播到所有 worker
    alpha_post_x = post_params_x.alpha_post;
    beta_post_x = post_params_x.beta_post;
    a_post_x = post_params_x.a_post;
    n_post_x = post_params_x.n_post;

    alpha_post_z = post_params_z.alpha_post;
    beta_post_z = post_params_z.beta_post;
    a_post_z = post_params_z.a_post;
    n_post_z = post_params_z.n_post;

    parfor i = 1:M
        % a) 从逆伽玛后验分布中抽样方差
        sigma2_x = invgamrnd(alpha_post_x, beta_post_x);
        sigma2_z = invgamrnd(alpha_post_z, beta_post_z);

        % b) 根据抽样到的方差，从正态后验分布中抽样均值
        mu_x = normrnd(a_post_x, sqrt(sigma2_x / (n_post_x+1)));
        mu_z = normrnd(a_post_z, sqrt(sigma2_z / (n_post_z+1)));
        
        s1 = sqrt(sigma2_x);
        s2 = sqrt(sigma2_z);
        % c) 使用抽样到的参数计算CEP
        try
            if s1 > 1e-9 && s2 > 1e-9
                R_bpt_samples(i) = calculate_cep_plugin(mu_x, mu_z, s1, s2);
            else
                R_bpt_samples(i) = NaN;
            end
        catch
            R_bpt_samples(i) = NaN;
        end
    end
    fprintf('  贝叶斯参数抽样完成。\n\n');

    R_bpt_samples = R_bpt_samples(~isnan(R_bpt_samples));

    % --- 3. 汇总结果 ---
    
    % 点估计 (公式 39)
    results.R_ba = mean(R_bpt_samples);
    
    % 置信区间和上界
    alpha = 1 - confidence_level;
    results.CI = quantile(R_bpt_samples, [alpha / 2, 1 - alpha / 2]);
    results.UB = quantile(R_bpt_samples, 1 - alpha);
end

% --- 辅助函数：计算后验分布参数 ---
function post_params = get_posterior_params(data_test, data_supp)
    n_test = length(data_test);
    n_supp = length(data_supp);
    
    % 1. 先验参数 (由补充样本决定)
    % 采用与标准贝叶斯理论一致的定义，这比规范中的公式更清晰
    if n_supp > 1
        a_prior = mean(data_supp);
        alpha_prior = 1;
        beta_prior = 0.5 / n_supp * sum((data_supp - a_prior).^2);
    else % 如果补充样本很少，使用弱信息先验
        a_prior = 0;
        alpha_prior = 0.001;
        beta_prior = 0.001;
    end
    
    % 2. 试验样本的统计量
    mean_test = mean(data_test);
    
    % 3. 后验参数 (标准更新公式)
    alpha_post = alpha_prior + (n_test + 1) / 2;
    a_post = (a_prior + n_test * mean_test) / (n_test + 1);
    
    sum_sq_err_test = sum((data_test - mean_test).^2) / n_test;
    
    beta_post = beta_prior + 0.5 * sum_sq_err_test * n_test;

    post_params.n_post = n_test;
    post_params.a_post = a_post;
    post_params.alpha_post = alpha_post;
    post_params.beta_post = beta_post;
end