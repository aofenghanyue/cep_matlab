function [CI, UB] = calculate_ci_bootstrap(R_hat, params, n, confidence_level, M)
% =========================================================================
% 使用参数自助法计算CEP的置信区间和置信上界
% =========================================================================
%
% 涉及流程: GJB 6289-2008, 第7.2节
%
%--------------------------------------------------------------------------
    mu1_hat = params.mu1;
    mu2_hat = params.mu2;
    s1_hat = params.s1;
    s2_hat = params.s2;
    
    R_boot_samples = zeros(M, 1);
    
    % 创建一个并行池以加速计算 (如果有并行计算工具箱)
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool;
    end
    
    % M次循环
    parfor i = 1:M
        % 1. 从估计的分布中生成新的样本
        u_new = normrnd(mu1_hat, s1_hat, [n, 1]);
        v_new = normrnd(mu2_hat, s2_hat, [n, 1]);
        
        % 2. 对新样本计算其统计参数
        mu1_boot = mean(u_new);
        mu2_boot = mean(v_new);
        s1_boot = std(u_new);
        s2_boot = std(v_new);
        
        % 3. 对新样本计算CEP值
        % 增加 try-catch 块防止因数值问题(如std为0)导致整个循环失败
        try
            if s1_boot > 1e-9 && s2_boot > 1e-9
                R_boot_samples(i) = calculate_cep_plugin(mu1_boot, mu2_boot, s1_boot, s2_boot);
            else
                R_boot_samples(i) = NaN;
            end
        catch
            R_boot_samples(i) = NaN;
        end
    end
    
    % 移除计算失败的样本
    R_boot_samples = R_boot_samples(~isnan(R_boot_samples));
    
    % 4. 计算分位数
    alpha = 1 - confidence_level;
    
    CI = R_hat.^2./quantile(R_boot_samples, [1 - alpha / 2, alpha / 2]);
    UB = R_hat.^2./quantile(R_boot_samples, alpha);
end