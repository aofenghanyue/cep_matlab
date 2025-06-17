function disp_cep_bayesian_results(results)
    fprintf('----------- [贝叶斯方法] 计算结果 -----------\n\n');
    fprintf('[CEP 贝叶斯点估计]\n');
    fprintf('  点估计 (R_ba): %.4f\n\n', results.R_ba);
    fprintf('[CEP 贝叶斯置信区间/上界]\n');
    fprintf('  置信区间: [%.4f, %.4f]\n', results.CI(1), results.CI(2));
    fprintf('  置信上界: %.4f\n\n', results.UB);
end