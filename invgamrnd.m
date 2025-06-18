function r = invgamrnd(alpha, beta, sz)
% =========================================================================
% 逆伽玛分布 (Inverse-Gamma Distribution) 随机数生成器
% =========================================================================
%
% 功能:
%   生成来自逆伽玛分布的随机数。其参数化与GJB 6289-2008和
%   大多数贝叶斯文献中的定义一致。(省略了国标中的归一化常数，因为matlab会自己处理)
%   PDF: p(x) ∝ (1/x)^(alpha+1) * exp(-beta/x)
%
% 输入:
%   alpha - 形状参数 (Shape parameter)
%   beta  - 尺度参数 (Scale parameter)
%   sz    - (可选) 输出随机数矩阵的尺寸，如 [m, n]
%
% 输出:
%   r     - 生成的随机数
%
% 实现原理:
%   如果 Y ~ Gamma(alpha, 1/beta)，那么 1/Y ~ Inverse-Gamma(alpha, beta)。
%   MATLAB的gamrnd(A,B)中，B是尺度参数。
%
%--------------------------------------------------------------------------
    if nargin < 3
        sz = [1, 1]; % 默认为单个值
    end

    % 从对应的伽玛分布 Gamma(alpha, 1/beta) 中抽样
    % gamrnd的第二个参数是伽玛分布的尺度参数，它等于逆伽玛尺度参数的倒数。
    gamma_samples = gamrnd(alpha, 1/beta, sz);
    
    % 取倒数即为逆伽玛分布的样本
    r = 1 ./ gamma_samples;
end