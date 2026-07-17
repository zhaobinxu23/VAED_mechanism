alpha_list = 0:4;
result_dir = 'batch_results_vaccine_variant';
if ~exist(result_dir, 'dir')
    mkdir(result_dir);
end
for ialpha = 1:length(alpha_list)
alpha = alpha_list(ialpha);
%% 1. 参数设置  second infection
p.N  = 400;      % 感染龄的分段数量 (bins)
p.da = 1.0;      % 感染龄步长 (Delta a, 单位通常是小时或步长单位)
p.dt_step = 1;   % 仅仅用于定义物理上的流速，通常为1(时间流逝速率)
p.infection_off = false;

% 感染终止阈值
p.extinction_threshold = 1;

%%
p.para(1) = 1e-20; % environmental antigen kon
p.para(2) = 0.5;% environmental antigen koff
p.para(3) = (2.2e17+1e13); % replenish constant pi 1 
p.para(4) = 0.5e13*0.5;% replenish constant pi 2
p.para(5) = 0.01; % decay constant of BCR IgM
p.para(6) = 0.005;% decay constant of BCR IgG
p.para(7) = 1; % k2  feedback constant of enviromental antigen-antibody complex
p.para(8) = 5e-7;% k2' feedback constant on PC cell regeneration
p.para(9) = 0.1; % decay constant of plasma Cell IgM
p.para(10) = 0.05;%% decay constant of plasma Cell IgG
p.para(11) = 4.4e9;% production constant of IgM
p.para(12) = 1.2e10;% production constant of IgG

p.para(13) = 0.05;% decay constant of IgM
p.para(14) = 0.025;% decay constant of IgG
p.para(15) = 5e4;% amplification constant of virus antigen
p.para(16) = 0.02;% transformation constant from IgM to IgG memory cell
p.para(17) = 0.1;% maximal production percentage of plasma cell 
p.para(18) = 0.5;% decay constant of complex  0.5
p.para(19) = 0; % virus replication constant
p.para(20) = 1e5;
p.para(21) = 1;

p.para_new(1) = 1e-22; 
p.para_new(2) = 1e-21;
p.para_new(3) = 1e-20; 
p.para_new(4) = 1e-19;
p.para_new(5) = 1e-18;
p.para_new(6) = 1e-17;
p.para_new(7) = 1e-16; 
p.para_new(8) = 1e-15;
p.para_new(9) = 1e-14; 
p.para_new(10) = 1e-13;

p.para_new_1(1) = 1e0; 
p.para_new_1(2) = 1e1;
p.para_new_1(3) = 1e2; 

p.para_new_1(4) = 1e3;
p.para_new_1(5) = 1e4; 
p.para_new_1(6) = 1e5;
p.para_new_1(7) = 1e6; 
p.para_new_1(8) = 1e7;
p.para_new_1(9) = 1e8; 
p.para_new_1(10) = 1e9;

% --- 病毒与易感细胞参数 ---
p.k4 = 1.0e-3;   % 病毒入侵常数
p.km = 1.0e6;    % 半饱和常数 2.0e6
p.k6 = 1.0e8;       % 易感细胞细胞再生率
p.k7 = 0.01;     % 易感细胞细胞自然死亡率
p.c_clear = 0.1; % 胞外病毒自然降解速率
p.Tc_generation = 5e-2;%% 抗原-抗体复合物刺激生成Tc细胞的速率 5e-5

% --- 胞内动力学参数 ---
p.v_start = 1;   % 刚感染时刻胞内病毒量
p.k5 = 0.1;     % 胞内病毒复制速率 (指数增长系数)
p.a_vec = (0:p.N-1)' * p.da; % 感染龄向量 [0, 1, 2, ..., 399]'

% 预计算胞内病毒量 (假设只跟感染龄有关)，这是一个 N x 1 向量
p.vin_vec = p.v_start .* exp(p.k5 .* p.a_vec);

% --- 裂解阈值参数 (Hill函数参数) ---
p.Tc_binding = 1e-5; % forward binding constant between Tc cell and infected cell
p.theta_lysis = 1e20; % 自然裂解阈值
p.theta_adcc  = 1e10;  % ADCC 结合分阈值
p.theta_tc = 1e8; % Tc 裂解的阈值
p.n_hill      = 2;    % Hill 系数 (控制平滑度/陡峭度)
p.k_lysis_max = 1;  % 超过阈值后的最大自然裂解速率
p.k_adcc_max  = 1; % 超过阈值后的最大ADCC裂解速率
p.k_tc_max = 1;% 超过阈值后的最大ADCC裂解速率

% --- Tc 细胞参数 (示例) ---
p.k_kill_tc = 0.1;   % Tc 杀伤系数

%% 2. 初始条件 (Initial Conditions)
% 状态向量结构 y = [T; I_1; ...; I_N; V; Tc; x0]
% 总长度 = 1 + N + 1 + 1 + 1401 = 1804

T0 = 1e10;                % 初始易感细胞
I0 = zeros(p.N, 1);      % 初始感染细胞 (各龄均为0)
V0 = 1e1;                % 初始胞外病毒
Tc0 = interp1(t_vac,y_vac(:,403),1000);               % 初始 Tc 细胞 1e5

%% 关于antibody部分：
x0 = zeros(3604,1);
%%  IgM distribution
mu = 9; % 均值
sigma = 0.80; % 标准差


% prob_A(1) = 0.001;
% prob_A(2:9) = 0;
prob_A(1) = normcdf(5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end
for i = 6:10
prob_A(i) = prob_A(11-i);
end

for i = 1:10
    for j = 1:10
        AA(10*(i-1)+j) = prob_A(i)*prob_A(j);
    end
end

M_1 = 1e15*0.5;
M_2 = 4e18*0.5;
G_1 = 1e15*0.5;
G_2 = 4e19*0.5;
E = 1e18;
C_1 = 1e13*0.5;
C_2 = 4e16*0.5;
C_3 = 1e13*0.5;
C_4 = 4e17*0.5;
P_M = 5e7*0.5;
P_G = 1e8*0.5;

%% *2
for i = 1:1800

     x0(i) = 1*interp1(t_vac,y_vac(:,i+403),1000);

end


for i = 1801:3600
     x0(i) = 10^(-alpha)*interp1(t_vac,y_vac(:,i+403),1000)+(1-10^(-alpha))*interp1(t_vac,y_vac(:,i+403),0);
end


 x0(3601) = E;
 x0(3602) = 1e1;% vaccine or virus
 x0(3603) = 0; % V-A1 Complex
 x0(3604) = 0; % V_A2 Complex

cum_death_0 = [0; 0; 0]; % 初始化这三个计数器为0


p.AA = AA;
p.AA_new = AA;
L_age_0 = zeros(p.N, 1);

y0 = [T0; ...
      I0; ...
      V0; ...
      Tc0; ...
      x0; ...
      cum_death_0; ...
      L_age_0];

%% 3. 运行模拟
t_span = [0, 1000]; % 模拟时间范围

% 使用 ode15s 求解，传入参数结构体 p
non_neg_indices = 1:length(y0);

% 修改 options
options = odeset('RelTol', 1e-4, ...
                 'AbsTol', 1e-8, ...
                 'NonNegative', non_neg_indices); % 强制所有变量非负

% 调用求解器

% [t1, y1] = ode15s(@(t,y) sys_ode_two_epitopes_new(t, y, p), t_span, y0, options);

[t_chal, y_chal, p_chal] = run_with_infection_extinction(y0, t_span, p, non_neg_indices);

challenge_file = fullfile(result_dir, ...
                sprintf('challenge_alpha_%d.mat', alpha));


 save(challenge_file, ...
                 't_chal', 'y_chal', 'p_chal');


end





function [t_all, y_all, p] = run_with_infection_extinction(y0, t_span, p, non_neg_indices)

    % 第一段：正常感染动力学，带 event 检测
    p.infection_off = false;

    options1 = odeset('RelTol', 1e-4, ...
                      'AbsTol', 1e-8, ...
                      'NonNegative', non_neg_indices, ...
                      'Events', @(t,y) infection_extinction_event(t, y, p));

    [t1, y1, te, ye, ie] = ode15s(@(t,y) sys_ode_two_epitopes_new_threshold(t, y, p), ...
                                  t_span, y0, options1);

    % 如果没有触发感染终止事件，则直接返回
    if isempty(te)
        t_all = t1;
        y_all = y1;
        return;
    end

    fprintf('Infection terminated at t = %.4f\n', te(end));

    % 触发感染终止后，清零相关变量
    y_reset = ye(end, :)';
    y_reset = apply_infection_extinction_reset(y_reset, p);

    % 第二段：感染终止后继续模拟剩余时间
    p.infection_off = true;

    t_span2 = [te(end), t_span(end)];

    options2 = odeset('RelTol', 1e-4, ...
                      'AbsTol', 1e-8, ...
                      'NonNegative', non_neg_indices);

    [t2, y2] = ode15s(@(t,y) sys_ode_two_epitopes_new_threshold(t, y, p), ...
                      t_span2, y_reset, options2);

    % 拼接两段结果，避免重复 te 那个时间点
    t_all = [t1; t2(2:end)];
    y_all = [y1; y2(2:end, :)];
end

function [value, isterminal, direction] = infection_extinction_event(t, y, p)

    N = p.N;

    I_total = sum(max(0, y(2:N+1)));
    V       = max(0, (y(N+2)+y(N+3606)));

    % 条件：V < 0.01 且 total infected cells < 0.01
    % 等价于 max(V, I_total) < 0.01
    value = max(V, I_total) - p.extinction_threshold;

    % 触发后终止当前积分
    isterminal = 1;

    % 只检测从正到负的穿越
    direction = -1;
end
function y = apply_infection_extinction_reset(y, p)

    N = p.N;

    % -----------------------------
    % 主感染模块
    % -----------------------------
    idx_I  = 2:N+1;
    idx_V  = N+2;

    y(idx_I) = 0;
    y(idx_V) = 0;

    % -----------------------------
    % 抗体模块 X
    % y(404:4007) 对应 X(1:3604)
    % -----------------------------
    X_start = N + 4;   % 对 N=400，就是 404

    % 第一套 epitope / antigen block


    % 第一套 virus-antibody complexes
    % X(801:1600) 对应 IgM-V, IgG-V 以及二级复合物相关变量
    idx_V_complex_1 = X_start + (801:1600) - 1;

    % 第二套 epitope / antigen block


    % 第二套 virus-antibody complexes
    idx_V_complex_2 = X_start + (2601:3400) - 1;

    % X(3602) 是抗体模块里面的 vaccine or virus antigen
    % X(3603), X(3604) 是 V-A1 complex 和 V-A2 complex
    idx_X_virus_and_complex = X_start + [3602, 3603, 3604] - 1;

    idx_zero = [idx_V_complex_1, ...
                idx_V_complex_2, ...
                idx_X_virus_and_complex];

    y(idx_zero) = 0;

    % 避免数值残留负数
    y = max(0, y);
end