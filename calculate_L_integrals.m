function L = calculate_L_integrals(type, powers, params, R)
% =========================================================================
% 计算 GJB 6289-2008 附录B 中的核心积分项 L1 和 L2
% =========================================================================
%
% 输入:
%   type   - 积分类型: 1 for L1, 2 for L2
%   powers - 幂次向量: 
%            for L1: [j, k] -> cos(phi)^j * sin(phi)^k
%            for L2: [i, j, k] -> r^i * cos(phi)^j * sin(phi)^k
%   params - 包含 mu1, mu2, s1, s2 的结构体
%   R      - CEP点估计值 R_hat (仅L1需要)
%
% 输出:
%   L      - 积分计算结果
%
%--------------------------------------------------------------------------
    
    mu1 = params.mu1;
    mu2 = params.mu2;
    s1 = params.s1;
    s2 = params.s2;

    % 定义附录B中的 a, b (注意这里的a,b与附录A不同)
    a = (1/(4*s2^2)) - (1/(4*s1^2));
    b = (1/(4*s2^2)) + (1/(4*s1^2));

    % 定义指数部分的核心表达式
    exp_term_fun = @(r, phi) exp(-b .* r.^2 + a .* r.^2 .* cos(2*phi) + ...
        (r/s1^2) .* mu1 .* cos(phi) + (r/s2^2) .* mu2 .* sin(phi));

    if type == 1 % 计算 L1(j, k)
        j = powers(1);
        k = powers(2);
        
        % L1的被积函数
        integrand_L1 = @(phi) (1/(2*pi)) .* (cos(phi).^j) .* (sin(phi).^k) .* exp_term_fun(R, phi);
        
        L = integral(integrand_L1, 0, 2*pi);
        
    elseif type == 2 % 计算 L2(i, j, k)
        i = powers(1);
        j = powers(2);
        k = powers(3);
        
        % L2的被积函数
        integrand_L2 = @(r, phi) (1/(2*pi)) .* (r.^i) .* (cos(phi).^j) .* (sin(phi).^k) .* exp_term_fun(r, phi);
        
        % 对 L2 进行双重积分
        L = integral2(integrand_L2, 0, R, 0, 2*pi);
        
    else
        error('未知的积分类型。请选择 1 (L1) 或 2 (L2)。');
    end
end