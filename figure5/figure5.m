
clear; clc; close all;

%% 1. 参数设置  first infection
p.N  = 400;      % 感染龄的分段数量 (bins)
p.da = 1.0;      % 感染龄步长 (Delta a, 单位通常是小时或步长单位)
p.dt_step = 1;   % 仅仅用于定义物理上的流速，通常为1(时间流逝速率)
% 是否已经感染终止
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
V0 = 10;                % 初始胞外病毒
Tc0 = 1e5;               % 初始 Tc 细胞 1e5

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
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j) = M_1*prob_A(i)*prob_A(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+100) = M_2*prob_A(i)*prob_A(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+200) = G_1*prob_A(i)*prob_A(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+300) = G_2*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+400) = C_1*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+500) = C_2*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+600) = C_3*prob_A(i)*prob_A(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+700) = C_4*prob_A(i)*prob_A(j);
    end
end

for i = 801:1600
   
 x0(i) = 0;
  
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1600) = P_M*prob_A(i)*prob_A(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1700) = P_G*prob_A(i)*prob_A(j);
    end
end

%%  precusor limitation model
mu = 9; % 均值
sigma = 0.80; % 标准差
prob_A_new(1) = normcdf(5, mu, sigma);
% 计算累积分布概率
for i = 2:5
prob_A_new(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end
for i = 6:10
prob_A_new(i) = prob_A_new(11-i);
end

for i = 1:10
    for j = 1:10
        AA_new(10*(i-1)+j) = prob_A_new(i)*prob_A_new(j);
    end
end

%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1800) = M_1*prob_A_new(i)*prob_A_new(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+100+1800) = M_2*prob_A_new(i)*prob_A_new(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+200+1800) = G_1*prob_A_new(i)*prob_A_new(j);
    end
end
%% *2
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+300+1800) = G_2*prob_A_new(i)*prob_A_new(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+400+1800) = C_1*prob_A_new(i)*prob_A_new(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+500+1800) = C_2*prob_A_new(i)*prob_A_new(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+600+1800) = C_3*prob_A_new(i)*prob_A_new(j);
    end
end
for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+700+1800) = C_4*prob_A_new(i)*prob_A_new(j);
    end
end

for i = 801+1800:1600+1800
   
 x0(i) = 0;
  
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1600+1800) = P_M*prob_A_new(i)*prob_A_new(j);
    end
end

for i = 1:10
    for j = 1:10
        x0(10*(i-1)+j+1700+1800) = P_G*prob_A_new(i)*prob_A_new(j);
    end
end

 x0(3601) = E;
 x0(3602) = 10;% vaccine or virus
 x0(3603) = 0; % V-A1 Complex
 x0(3604) = 0; % V_A2 Complex

cum_death_0 = [0; 0; 0]; % 初始化这三个计数器为0


p.AA = AA;
p.AA_new = AA_new;

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
[t, y, p] = run_with_infection_extinction(y0, t_span, p, non_neg_indices);

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
Tc0 = interp1(t,y(:,403),0);               % 初始 Tc 细胞 1e5

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

     x0(i) = 1*interp1(t,y(:,i+403),1000);

end


for i = 1801:3600
     x0(i) = 1*interp1(t,y(:,i+403),0);
end


 x0(3601) = E;
 x0(3602) = 1e1;% vaccine or virus
 x0(3603) = 0; % V-A1 Complex
 x0(3604) = 0; % V_A2 Complex

cum_death_0 = [0; 0; 0]; % 初始化这三个计数器为0


p.AA = AA;
p.AA_new = AA_new;
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

[t1, y1, p] = run_with_infection_extinction(y0, t_span, p, non_neg_indices);




%% ============================================================
%  Figure 4a-f plotting code only
%  Uniform absolute time grid: 0:1:1000
%  Keep main simulation and ODE unchanged
%  ============================================================

%% ------------------------------------------------------------
%  0. Uniform time grid and basic indices
%  ------------------------------------------------------------

t_grid = (0:1:1000)';

N = p.N;
age_vec = p.a_vec(:);

idx_T  = 1;
idx_I  = 2:(N+1);
idx_V  = N + 2;
idx_Tc = N + 3;

idx_x0_start = N + 4;
idx_x0_end   = idx_x0_start + 3604 - 1;

idx_X3603 = idx_x0_start + 3603 - 1;
idx_X3604 = idx_x0_start + 3604 - 1;

idx_cum_nat  = idx_x0_end + 1;
idx_cum_adcc = idx_x0_end + 2;
idx_cum_tc   = idx_x0_end + 3;

idx_L_age_start = idx_x0_end + 3 + 1;
idx_L_age_end   = idx_L_age_start + N - 1;

%% ------------------------------------------------------------
%  1. Interpolate full solutions to the same absolute time grid
%  ------------------------------------------------------------

% interp1 对矩阵 y 会逐列插值，所以每个状态变量都会被插值到 t_grid
Y_first_grid = interp1(t, y, t_grid, 'linear', 'extrap');
Y_second_grid = interp1(t1, y1, t_grid, 'linear', 'extrap');

Y_first_grid = real(Y_first_grid);
Y_second_grid = real(Y_second_grid);

Y_first_grid(~isfinite(Y_first_grid)) = 0;
Y_second_grid(~isfinite(Y_second_grid)) = 0;

% 数值误差可能产生极小负数，统一截断
Y_first_grid(Y_first_grid < 0) = 0;
Y_second_grid(Y_second_grid < 0) = 0;

%% ------------------------------------------------------------
%  2. Extract variables on t_grid
%  ------------------------------------------------------------

V_first  = Y_first_grid(:, idx_V);
V_second = Y_second_grid(:, idx_V);

% infectious virus = free virus + X(3603)
infectious_V_first  = Y_first_grid(:, idx_V)  + Y_first_grid(:, idx_X3603);
infectious_V_second = Y_second_grid(:, idx_V) + Y_second_grid(:, idx_X3603);

I_first  = Y_first_grid(:, idx_I);
I_second = Y_second_grid(:, idx_I);

total_I_first  = sum(I_first, 2);
total_I_second = sum(I_second, 2);

cum_nat_first  = Y_first_grid(:, idx_cum_nat);
cum_adcc_first = Y_first_grid(:, idx_cum_adcc);
cum_tc_first   = Y_first_grid(:, idx_cum_tc);

cum_nat_second  = Y_second_grid(:, idx_cum_nat);
cum_adcc_second = Y_second_grid(:, idx_cum_adcc);
cum_tc_second   = Y_second_grid(:, idx_cum_tc);

cum_total_first  = cum_nat_first  + cum_adcc_first  + cum_tc_first;
cum_total_second = cum_nat_second + cum_adcc_second + cum_tc_second;

%% ============================================================
%  Figure 4a
%  Free extracellular virus and infectious virus
%  ============================================================

figure('Color','w','Position',[100 100 1150 480]);

subplot(1,2,1);
hold on;
plot(t_grid, V_first, 'k-', 'LineWidth', 2);
plot(t_grid, V_second, 'r-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Free extracellular virus');
title('Free extracellular virus');
legend({'Primary infection','Secondary variant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

subplot(1,2,2);
hold on;
plot(t_grid, infectious_V_first, 'k-', 'LineWidth', 2);
plot(t_grid, infectious_V_second, 'r-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Infectious virus, V + X_{3603}');
title('Infectious virus');
legend({'Primary infection','Secondary variant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

sgtitle('a, Extracellular viral load during primary and secondary infection');

saveas(gcf, 'figure4a_viral_load_free_and_infectious.png');
saveas(gcf, 'figure4a_viral_load_free_and_infectious.fig');

%% ============================================================
%  Figure 4b
%  Total infected-cell burden over time
%  ============================================================

figure('Color','w','Position',[100 100 650 500]);
hold on;
plot(t_grid, total_I_first, 'k-', 'LineWidth', 2);
plot(t_grid, total_I_second, 'r-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Total infected cells');
title('b, Total infected-cell burden over time');
legend({'Primary infection','Secondary variant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

saveas(gcf, 'figure4b_total_infected_cell_burden.png');
saveas(gcf, 'figure4b_total_infected_cell_burden.fig');

%% ============================================================
%  Figure 4c
%  Cumulative infected-cell lysis
%  ============================================================

figure('Color','w','Position',[100 100 1200 520]);

subplot(1,2,1);
hold on;
plot(t_grid, cum_nat_first,   'b-',  'LineWidth', 2);
plot(t_grid, cum_adcc_first,  'r-',  'LineWidth', 2);
plot(t_grid, cum_tc_first,    'g-',  'LineWidth', 2);
plot(t_grid, cum_total_first, 'k--', 'LineWidth', 2);
xlabel('Time after primary infection');
ylabel('Cumulative lysed infected cells');
title('Primary infection');
legend({'Natural lysis','ADCC lysis','Tc killing','Total'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

subplot(1,2,2);
hold on;
plot(t_grid, cum_nat_second,   'b-',  'LineWidth', 2);
plot(t_grid, cum_adcc_second,  'r-',  'LineWidth', 2);
plot(t_grid, cum_tc_second,    'g-',  'LineWidth', 2);
plot(t_grid, cum_total_second, 'k--', 'LineWidth', 2);
xlabel('Time after secondary infection');
ylabel('Cumulative lysed infected cells');
title('Secondary variant infection');
legend({'Natural lysis','ADCC lysis','Tc killing','Total'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

sgtitle('c, Cumulative infected-cell lysis');

saveas(gcf, 'figure4c_cumulative_infected_cell_lysis.png');
saveas(gcf, 'figure4c_cumulative_infected_cell_lysis.fig');

%% Total cumulative lysis comparison
figure('Color','w','Position',[100 100 650 500]);
hold on;
plot(t_grid, cum_total_first, 'k-', 'LineWidth', 2);
plot(t_grid, cum_total_second, 'r-', 'LineWidth', 2);
xlabel('Time after infection');
ylabel('Total cumulative lysed infected cells');
title('Total cumulative infected-cell lysis');
legend({'Primary infection','Secondary variant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

saveas(gcf, 'figure4c_total_cumulative_lysis_comparison.png');
saveas(gcf, 'figure4c_total_cumulative_lysis_comparison.fig');

%% ============================================================
%  Figure 4d
%  Distribution of infected-cell by infection age
%  ============================================================

Z_I_first  = log10(I_first' + 1);
Z_I_second = log10(I_second' + 1);

Z_I_first(~isfinite(Z_I_first)) = 0;
Z_I_second(~isfinite(Z_I_second)) = 0;

common_I_min = 0;
common_I_max = max([Z_I_first(:); Z_I_second(:)]);
if ~isfinite(common_I_max) || common_I_max <= common_I_min
    common_I_max = 1;
end

figure('Color','w','Position',[100 100 1200 500]);

subplot(1,2,1);
imagesc(t_grid, age_vec, Z_I_first);
set(gca, 'YDir', 'normal');
caxis([common_I_min common_I_max]);
colormap hot;
cb1 = colorbar;
ylabel(cb1, 'log_{10}(infected cells + 1)');
xlabel('Time after primary infection');
ylabel('Infection age');
title('Primary infection');
set(gca, 'FontSize', 12);

subplot(1,2,2);
imagesc(t_grid, age_vec, Z_I_second);
set(gca, 'YDir', 'normal');
caxis([common_I_min common_I_max]);
colormap hot;
cb2 = colorbar;
ylabel(cb2, 'log_{10}(infected cells + 1)');
xlabel('Time after secondary infection');
ylabel('Infection age');
title('Secondary variant infection');
set(gca, 'FontSize', 12);

sgtitle('d, Distribution of infected cells by infection age');

saveas(gcf, 'figure4d_infected_cell_distribution_by_age.png');
saveas(gcf, 'figure4d_infected_cell_distribution_by_age.fig');

%% ============================================================
%  Figure 4e and 4f
%  Lysis distribution and viral release on absolute time grid
%  ============================================================

if size(Y_first_grid,2) >= idx_L_age_end && size(Y_second_grid,2) >= idx_L_age_end

    %% ------------------------------------------------------------
    %  4e. Extract cumulative L_age on t_grid
    %  ------------------------------------------------------------

    L_age_first_cum  = Y_first_grid(:, idx_L_age_start:idx_L_age_end);
    L_age_second_cum = Y_second_grid(:, idx_L_age_start:idx_L_age_end);

    % 在统一时间轴上求差分，因为 t_grid 间隔为 1
    % L_age_rate(k,:) 表示 t_grid(k-1) 到 t_grid(k) 之间的裂解数
    L_age_first_rate = [zeros(1, N); diff(L_age_first_cum, 1, 1)];
    L_age_second_rate = [zeros(1, N); diff(L_age_second_cum, 1, 1)];

    L_age_first_rate = real(L_age_first_rate);
    L_age_second_rate = real(L_age_second_rate);

    L_age_first_rate(~isfinite(L_age_first_rate)) = 0;
    L_age_second_rate(~isfinite(L_age_second_rate)) = 0;

    L_age_first_rate(L_age_first_rate < 0) = 0;
    L_age_second_rate(L_age_second_rate < 0) = 0;

    Z_L_first  = log10(L_age_first_rate' + 1);
    Z_L_second = log10(L_age_second_rate' + 1);

    Z_L_first(~isfinite(Z_L_first)) = 0;
    Z_L_second(~isfinite(Z_L_second)) = 0;

    common_L_min = 0;
    common_L_max = max([Z_L_first(:); Z_L_second(:)]);
    if ~isfinite(common_L_max) || common_L_max <= common_L_min
        common_L_max = 1;
    end

    figure('Color','w','Position',[100 100 1200 500]);

    subplot(1,2,1);
    imagesc(t_grid, age_vec, Z_L_first);
    set(gca, 'YDir', 'normal');
    caxis([common_L_min common_L_max]);
    colormap hot;
    cb1 = colorbar;
    ylabel(cb1, 'log_{10}(lysed cells per time step + 1)');
    xlabel('Time after primary infection');
    ylabel('Infection age');
    title('Primary infection');
    set(gca, 'FontSize', 12);

    subplot(1,2,2);
    imagesc(t_grid, age_vec, Z_L_second);
    set(gca, 'YDir', 'normal');
    caxis([common_L_min common_L_max]);
    colormap hot;
    cb2 = colorbar;
    ylabel(cb2, 'log_{10}(lysed cells per time step + 1)');
    xlabel('Time after secondary infection');
    ylabel('Infection age');
    title('Secondary variant infection');
    set(gca, 'FontSize', 12);

    sgtitle('e, Distribution of infected-cell lysis by infection age');

    saveas(gcf, 'figure4e_lysis_distribution_by_infection_age.png');
    saveas(gcf, 'figure4e_lysis_distribution_by_infection_age.fig');

    %% Mean infection age at lysis
    age_row = age_vec(:)';

    mean_lysis_age_first = sum(L_age_first_rate .* age_row, 2) ./ ...
                           max(sum(L_age_first_rate, 2), eps);

    mean_lysis_age_second = sum(L_age_second_rate .* age_row, 2) ./ ...
                            max(sum(L_age_second_rate, 2), eps);

    mean_lysis_age_first(~isfinite(mean_lysis_age_first)) = 0;
    mean_lysis_age_second(~isfinite(mean_lysis_age_second)) = 0;

    figure('Color','w','Position',[100 100 650 500]);
    hold on;
    plot(t_grid, mean_lysis_age_first, 'k-', 'LineWidth', 2);
    plot(t_grid, mean_lysis_age_second, 'r-', 'LineWidth', 2);
    xlabel('Time after infection');
    ylabel('Mean infection age at lysis');
    title('Mean infection age of lysed cells');
    legend({'Primary infection','Secondary variant infection'}, 'Location','best');
    box on;
    set(gca, 'FontSize', 12);

    saveas(gcf, 'figure4e_mean_lysis_age_comparison.png');
    saveas(gcf, 'figure4e_mean_lysis_age_comparison.fig');

    %% ------------------------------------------------------------
    %  4f. Viral release per lysed cell
    %  ------------------------------------------------------------

    vin_row = p.vin_vec(:)';

    viral_release_first_age = L_age_first_rate .* vin_row;
    viral_release_second_age = L_age_second_rate .* vin_row;

    viral_release_first_age(~isfinite(viral_release_first_age)) = 0;
    viral_release_second_age(~isfinite(viral_release_second_age)) = 0;

    Z_release_first  = log10(viral_release_first_age' + 1);
    Z_release_second = log10(viral_release_second_age' + 1);

    Z_release_first(~isfinite(Z_release_first)) = 0;
    Z_release_second(~isfinite(Z_release_second)) = 0;

    common_R_min = 0;
    common_R_max = max([Z_release_first(:); Z_release_second(:)]);
    if ~isfinite(common_R_max) || common_R_max <= common_R_min
        common_R_max = 1;
    end

    figure('Color','w','Position',[100 100 1200 500]);

    subplot(1,2,1);
    imagesc(t_grid, age_vec, Z_release_first);
    set(gca, 'YDir', 'normal');
    caxis([common_R_min common_R_max]);
    colormap hot;
    cb1 = colorbar;
    ylabel(cb1, 'log_{10}(virions released per time step + 1)');
    xlabel('Time after primary infection');
    ylabel('Infection age');
    title('Primary infection');
    set(gca, 'FontSize', 12);

    subplot(1,2,2);
    imagesc(t_grid, age_vec, Z_release_second);
    set(gca, 'YDir', 'normal');
    caxis([common_R_min common_R_max]);
    colormap hot;
    cb2 = colorbar;
    ylabel(cb2, 'log_{10}(virions released per time step + 1)');
    xlabel('Time after secondary infection');
    ylabel('Infection age');
    title('Secondary variant infection');
    set(gca, 'FontSize', 12);

    sgtitle('f, Viral release by lysed infected cells');

    saveas(gcf, 'figure4f_viral_release_by_infection_age.png');
    saveas(gcf, 'figure4f_viral_release_by_infection_age.fig');

    %% Viral release per lysed cell, weighted average over infection age
    release_per_lysed_first = sum(viral_release_first_age, 2) ./ ...
                              max(sum(L_age_first_rate, 2), eps);

    release_per_lysed_second = sum(viral_release_second_age, 2) ./ ...
                               max(sum(L_age_second_rate, 2), eps);

    release_per_lysed_first(~isfinite(release_per_lysed_first)) = 0;
    release_per_lysed_second(~isfinite(release_per_lysed_second)) = 0;

    figure('Color','w','Position',[100 100 650 500]);
    hold on;
    plot(t_grid, release_per_lysed_first, 'k-', 'LineWidth', 2);
    plot(t_grid, release_per_lysed_second, 'r-', 'LineWidth', 2);
    set(gca, 'YScale', 'log');
    xlabel('Time after infection');
    ylabel('Viral release per lysed cell');
    title('Viral release per lysis event');
    legend({'Primary infection','Secondary variant infection'}, 'Location','best');
    box on;
    set(gca, 'FontSize', 12);

    saveas(gcf, 'figure4f_viral_release_per_lysed_cell.png');
    saveas(gcf, 'figure4f_viral_release_per_lysed_cell.fig');

else
    warning('L_age variables are not found in y/y_second. Figure 4e and 4f are skipped.');
end

%% ============================================================
%  Save plotted data on uniform time grid
%  ============================================================

save('figure4abcdef_plot_data_uniform_time_grid.mat', ...
     't_grid', ...
     'Y_first_grid', 'Y_second_grid', ...
     'V_first', 'V_second', ...
     'infectious_V_first', 'infectious_V_second', ...
     'total_I_first', 'total_I_second', ...
     'cum_nat_first', 'cum_adcc_first', 'cum_tc_first', ...
     'cum_nat_second', 'cum_adcc_second', 'cum_tc_second', ...
     'cum_total_first', 'cum_total_second');



%%
primarySim.t = t;
primarySim.y = y;
primarySim.p = p;

elisaTimeGrid = 0:10:max(t);

primaryElisa = compute_elisa_signal_two_epitopes_primary( ...
    primarySim, ...
    elisaTimeGrid, ...
    'Primary infection');

save('primary_infection_elisa_neutralizing_vs_nonneutralizing.mat', ...
    'primaryElisa', ...
    'elisaTimeGrid');


primarySim.t = t1;
primarySim.y = y1;
primarySim.p = p;

elisaTimeGrid = 0:10:max(t1);

SecondaryElisa = compute_elisa_signal_two_epitopes_primary( ...
    primarySim, ...
    elisaTimeGrid, ...
    'Primary infection');

save('Secondary_infection_elisa_neutralizing_vs_nonneutralizing.mat', ...
    'SecondaryElisa', ...
    'elisaTimeGrid');

 plot_primary_infection_elisa_neu_vs_nonNeu_new( ...
    primaryElisa, ...
    SecondaryElisa, ...
    'primary_vs_secondary_elisa_neutralizing_vs_nonneutralizing');
%

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