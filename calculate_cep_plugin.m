function R = calculate_cep_plugin(varargin)
% =========================================================================
% 计算CEP的代入型点估计值 R_hat
% --- 在这里选择您想使用的计算方法 ---
% 1: 笛卡尔坐标积分 (不推荐)
% 2: 极坐标积分 (稳定, 推荐)
method = 1;
if method ==1
    R = calculate_cep_plugin_circle(varargin{:});
else
    R = calculate_cep_plugin_polar(varargin{:});
end
end