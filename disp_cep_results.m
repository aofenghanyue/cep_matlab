function disp_cep_results(results, method_name)
    fprintf('----------- [%s] 计算结果 -----------\n\n', method_name);
    fprintf('[基础统计参数]\n');
    fprintf('  原始数据 (变换前):\n');
    fprintf('    μ_x = %.4f, σ_x = %.4f\n', results.basic_stats.mux, results.basic_stats.sx);
    fprintf('    μ_z = %.4f, σ_z = %.4f\n', results.basic_stats.muz, results.basic_stats.sz);
    fprintf('    ρ   = %.5f\n', results.basic_stats.rho);
    if ~results.tests.independence.is_independent
        fprintf('  由于x与z独立性检验未通过, 已进行正交变换。\n');
        fprintf('  变换后用于计算的参数 (独立):\n');
        fprintf('    μ₁ = %.4f, σ₁ = %.4f\n', results.params.mu1, results.params.s1);
        fprintf('    μ₂ = %.4f, σ₂ = %.4f\n', results.params.mu2, results.params.s2);
    end
    fprintf('\n');
    fprintf('[CEP 点估计 (Point Estimation)]\n');
    fprintf('  代入型点估计 (Plug-in): R̂ = %.4f\n\n', results.point_estimate.R_hat);
    fprintf('[CEP 区间估计 (Interval Estimation)]\n');
    fprintf('  方法 1: 参数自助法 (Parametric Bootstrap) - 推荐\n');
    fprintf('    置信区间: [%.4f, %.4f]\n', results.ci_bootstrap.CI(1), results.ci_bootstrap.CI(2));
    fprintf('    置信上界: %.4f\n\n', results.ci_bootstrap.UB);
    fprintf('  方法 2: 一阶逼近法 (First-Order Approximation)\n');
    fprintf('    置信区间: [%.4f, %.4f]\n', results.ci_first_order.CI(1), results.ci_first_order.CI(2));
    fprintf('    置信上界: %.4f\n\n', results.ci_first_order.UB);
end