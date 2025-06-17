function R = calculate_cep_plugin_polar(varargin)
% =========================================================================
% 计算CEP的代入型点估计值 R_hat
% =========================================================================
%   GJB 6289-2008 附录A.2的极坐标积分方法。
%   避免了笛卡尔坐标积分的奇异性问题，提高了计算的稳定性和效率。
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

    % --- 1. 定义极坐标下的被积函数 (不含雅可比行列式r) ---
    % integrand_polar(r, phi)
    pdf_polar = @(r, phi) (1 / (2 * pi * s1 * s2)) * ...
        exp(-0.5 * ( ((r .* cos(phi) - mu1).^2 / s1^2) + ...
                      ((r .* sin(phi) - mu2).^2 / s2^2) ));

    % --- 2. 定义内层积分 (对角度 phi 积分) ---
    % 该函数计算在半径 r 处的概率密度环的值乘以 r (因为 dr dphi -> r dr dphi)
    inner_integral = @(r) r .* integral(@(phi) pdf_polar(r, phi), 0, 2*pi, 'ArrayValued', true);
    
    % --- 3. 定义需要求根的函数 f(R) = integral(inner_integral) - 0.5 ---
    % 外层积分，从 0 到 R
    fun_to_zero = @(R) integral(inner_integral, 0, R) - 0.5;

    % --- 4. 求解 ---
    % 提供一个合理的初始猜测值
    initial_guess = sqrt(s1^2 + s2^2 + mu1^2 + mu2^2);
    if initial_guess < 1e-6
        initial_guess = 1;
    end
    
    % 使用 fzero 求解 R
    options = optimset('Display', 'off');
    try
        R = fzero(fun_to_zero, initial_guess, options);
    catch ME
        % 如果 fzero 失败 (例如初始点两侧函数值同号)，尝试更大的搜索范围
        if strcmp(ME.identifier, 'MATLAB:fzero:ValuesAtEndPtsSameSign')
            warning('fzero在初始点附近失败，尝试扩大搜索范围。');
            % 创建一个更宽的搜索区间
            search_interval = [0, initial_guess * 5]; 
            R = fzero(fun_to_zero, search_interval, options);
        else
            rethrow(ME);
        end
    end
end