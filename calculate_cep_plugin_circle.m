function R = calculate_cep_plugin_circle(varargin)
% =========================================================================
% 计算CEP的代入型点估计值 R_hat
% =========================================================================
%
% 功能:
%   根据给定的统计参数，求解CEP积分方程，得到CEP值。
%
% 输入:
%   % 允许多种输入方式
%   1. 输入一个结构体: calculate_cep_plugin(params)
%   2. 输入四个独立参数: calculate_cep_plugin(mu1, mu2, s1, s2)
%
% 输出:
%   R      - CEP的点估计值
%
% 涉及公式: GJB 6289-2008, 公式(11)
%
%--------------------------------------------------------------------------
    if nargin == 1
        params = varargin{1};
        mu1 = params.mu1;
        mu2 = params.mu2;
        s1 = params.s1;
        s2 = params.s2;
    elseif nargin == 4
        mu1 = varargin{1};
        mu2 = varargin{2};
        s1 = varargin{3};
        s2 = varargin{4};
    else
        error('输入参数数量不正确。请输入一个params结构体或四个独立的参数(mu1, mu2, s1, s2)。');
    end

    % 定义二维正态分布的概率密度函数 (PDF)
    pdf = @(x, y) (1 / (2 * pi * s1 * s2)) * ...
                   exp(-0.5 * (((x - mu1).^2 / s1^2) + ((y - mu2).^2 / s2^2)));
    
    % 定义需要求根的函数 f(R) = integral(PDF) - 0.5
    fun_to_zero = @(R) integral2(pdf, -R, R, @(x) -sqrt(R^2 - x.^2), @(x) sqrt(R^2 - x.^2)) - 0.5;
    
    % 提供一个合理的初始猜测值
    initial_guess = sqrt(s1^2 + s2^2 + mu1^2 + mu2^2);
    if initial_guess < 1e-6
        initial_guess = 1;
    end
    
    % 使用 fzero 求解 R
    options = optimset('Display', 'off');
    R = fzero(fun_to_zero, initial_guess, options);
end