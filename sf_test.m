function [h, p] = sf_test(x, alpha)
% =========================================================================
% Shapiro-Francia 正态性检验 (sf_test)
% =========================================================================
%
% 功能:
%   检验输入数据向量x是否来自正态分布。这是对Shapiro-Wilk检验的一个
%   有效且易于实现的替代方案，不依赖于任何工具箱。
%   计算步骤：
%   1.对输入的样本数据 x 进行升序排序，得到 x_sorted。
%   2.计算这组数据对应的"正态得分"（理论上的标准正态分位数）。一个常用的近似计算方法是Blom变换：m = norminv((i - 3/8) / (n + 1/4))，其中 i 是排序后的序号，n 是样本量。
%   3.计算 x_sorted 和 m 之间的皮尔逊相关系数 corr(x_sorted, m)。
%   4.检验统计量 W' 就是这个相关系数的平方：W' = (corr(x_sorted, m))^2。
%   5.W' 的值越接近1，说明数据越接近正态分布。
%   6.为了得到p值，我们需要将 W' 转换为一个近似服从标准正态分布的Z分数。Royston (1983) 提出了一个有效的转换公式。
%   7.根据 Z 分数，计算出p值。
%
% 输入:
%   x     - 待检验的数据向量
%   alpha - 显著性水平 (可选, 默认为 0.05)
%
% 输出:
%   h     - 检验结果: 0表示不能拒绝原假设(正态), 1表示拒绝(非正态)
%   p     - 检验的P值
%
% 参考文献:
%   - Royston, J. P. (1983). "A simple method for evaluating the
%     Shapiro-Francia W' test of non-normality." Journal of the Royal
%     Statistical Society: Series D (The Statistician), 32(3), 297-300.
%
%--------------------------------------------------------------------------

    if nargin < 2, alpha = 0.10; end
    x = x(:); % 将输入的 x (无论是行是列) 统一转换成列向量
    x = x(~isnan(x)); % 移除NaN值
    x = sort(x);      % 升序排序
    n = length(x);

    if n < 4 || n > 5000
        warning('Shapiro-Francia检验适用于样本量在4到5000之间。');
    end

    % 1. 计算正态得分 (Normal Scores) - Blom's plotting position
    i = (1:n)';
    m = norminv((i - 3/8) / (n + 1/4));

    % 2. 计算检验统计量 W'
    % W' 是数据与正态得分之间相关系数的平方
    W_prime = (corr(x, m))^2;

    % 3. 将 W' 转换为近似正态分布的 Z 分数 (Royston's transform)
    nu = log(n);
    u = log(nu) - nu;
    mu = -1.2725 + 1.0521 * u;
    sigma = 1.0308 - 0.26758 * u;
    
    g = log(1 - W_prime);
    z = (g - mu) / sigma;

    % 4. 计算 P 值
    % 这是双侧检验，但我们通常关心的是 "不够正态" 的一侧，
    % 即 W' 过小，g 过大，z 过大。所以我们计算右尾概率。
    p = 1 - normcdf(z);

    % 5. 确定检验结果 h
    if p < alpha
        h = 1; % 拒绝原假设，数据非正态
    else
        h = 0; % 不能拒绝原假设，数据可视为正态
    end
end