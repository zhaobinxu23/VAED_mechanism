clear; clc; close all;

%% 1. 参数设置  first infection
p.N  = 400;      % 感染龄的分段数量 (bins)
p.da = 1.0;      % 感染龄步长 (Delta a, 单位通常是小时或步长单位)
p.dt_step = 1;   % 仅仅用于定义物理上的流速，通常为1(时间流逝速率)

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
p.para(15) = 1e4;% amplification constant of virus antigen
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
p.km = 2.0e6;    % 半饱和常数 2.0e6
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
p.theta_lysis = 1e12; % 自然裂解阈值
p.theta_adcc  = 1e8;  % ADCC 结合分阈值
p.theta_tc = 1e7; % Tc 裂解的阈值
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
sigma = 0.8; % 标准差


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
sigma = 0.8; % 标准差
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
y0 = [T0; I0; V0; Tc0; x0; cum_death_0]; 

%% 3. 运行模拟
t_span = [0, 1000]; % 模拟时间范围

% 使用 ode15s 求解，传入参数结构体 p
non_neg_indices = 1:length(y0);

% 修改 options
options = odeset('RelTol', 1e-4, ...
                 'AbsTol', 1e-8, ...
                 'NonNegative', non_neg_indices); % 强制所有变量非负



% 调用求解器
[t, y] = ode15s(@(t,y) sys_ode_two_epitopes_new(t, y, p), t_span, y0, options);
% Primary infection: neutralizing vs non-neutralizing antibody landscapes
% ========================================================================

%% ========================================================================
% Primary infection: log-scaled time-dependent heatmaps
% Neutralizing vs non-neutralizing antibodies
% ========================================================================

% -------------------------------------------------------------------------
% Time grid and affinity distribution
% -------------------------------------------------------------------------
t_grid = 0:2:1000;     % 如果你的模拟终点是 1000
% t_grid = 0:2:2000;   % 如果你的模拟终点是 2000，则使用这一行

dist2 = compute_affinity_distribution_two_epitopes(t, y, t_grid, p);

% -------------------------------------------------------------------------
% Use log-transformed abundance for visualization
% -------------------------------------------------------------------------
Neu_abs    = dist2.Neu_total_kd';       % time x Kd, absolute abundance
NonNeu_abs = dist2.NonNeu_total_kd';    % time x Kd, absolute abundance

Neu_plot    = log10(Neu_abs + 1);
NonNeu_plot = log10(NonNeu_abs + 1);

% Shared log color scale for direct comparison
clim_shared = [0, max([Neu_plot(:); NonNeu_plot(:)])];
if clim_shared(2) <= clim_shared(1)
    clim_shared = [0 1];
end

% -------------------------------------------------------------------------
% High-contrast colormap
% -------------------------------------------------------------------------
sci_cmap = interp1([1 64 128 192 256], ...
    [0.05 0.10 0.35; ...   % deep navy
     0.10 0.35 0.75; ...   % blue
     0.95 0.95 0.55; ...   % pale yellow
     0.95 0.45 0.15; ...   % orange
     0.55 0.00 0.00], ...  % dark red
    1:256);

set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesLineWidth', 1.0);

fig = figure('Color', 'w', ...
             'Units', 'centimeters', ...
             'Position', [2 2 20 10], ...
             'Renderer', 'painters');

tl = tiledlayout(fig, 1, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% ------------------------------------------------------------------------
% Panel a: RBD-specific neutralizing antibodies
% -------------------------------------------------------------------------
ax1 = nexttile(tl, 1);

imagesc(ax1, dist2.kd_values, t_grid, Neu_plot);
axis(ax1, 'xy');
axis(ax1, 'tight');

colormap(ax1, sci_cmap);
caxis(ax1, clim_shared);

set(ax1, ...
    'FontName', 'Arial', ...
    'FontSize', 8.5, ...
    'LineWidth', 0.9, ...
    'TickDir', 'out', ...
    'Layer', 'top', ...
    'Box', 'off', ...
    'Color', 'w', ...
    'XColor', [0.15 0.15 0.15], ...
    'YColor', [0.15 0.15 0.15]);

xlabel(ax1, 'Affinity group, log_{10}(K_D)', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold');

ylabel(ax1, 'Time after primary infection', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold');

title(ax1, {'RBD-specific', 'neutralizing antibodies'}, ...
    'FontSize', 10.5, ...
    'FontWeight', 'bold');

xticks(ax1, -31:3:-13);
yticks(ax1, 0:250:max(t_grid));

text(ax1, -0.16, 1.04, 'a', ...
    'Units', 'normalized', ...
    'FontSize', 13, ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% Panel b: non-RBD-specific non-neutralizing antibodies
% -------------------------------------------------------------------------
ax2 = nexttile(tl, 2);

imagesc(ax2, dist2.kd_values, t_grid, NonNeu_plot);
axis(ax2, 'xy');
axis(ax2, 'tight');

colormap(ax2, sci_cmap);
caxis(ax2, clim_shared);

set(ax2, ...
    'FontName', 'Arial', ...
    'FontSize', 8.5, ...
    'LineWidth', 0.9, ...
    'TickDir', 'out', ...
    'Layer', 'top', ...
    'Box', 'off', ...
    'Color', 'w', ...
    'XColor', [0.15 0.15 0.15], ...
    'YColor', [0.15 0.15 0.15]);

xlabel(ax2, 'Affinity group, log_{10}(K_D)', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold');

ylabel(ax2, 'Time after primary infection', ...
    'FontSize', 9.5, ...
    'FontWeight', 'bold');

title(ax2, {'non-RBD-specific', 'non-neutralizing antibodies'}, ...
    'FontSize', 10.5, ...
    'FontWeight', 'bold');

xticks(ax2, -31:3:-13);
yticks(ax2, 0:250:max(t_grid));

text(ax2, -0.16, 1.04, 'b', ...
    'Units', 'normalized', ...
    'FontSize', 13, ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------------------
% Shared colorbar
% -------------------------------------------------------------------------
cb = colorbar(ax2);
cb.Layout.Tile = 'east';
cb.Label.String = 'log_{10}(antibody abundance + 1)';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 9.5;
cb.Label.FontWeight = 'bold';
cb.FontName = 'Arial';
cb.FontSize = 8.5;
cb.LineWidth = 0.8;
cb.Box = 'off';

%% ------------------------------------------------------------------------
% Main title
% -------------------------------------------------------------------------
sgtitle(tl, ...
    {'Primary infection with ancestral virus', ...
     'RBD and non-RBD epitopes induce comparable antibody expansion'}, ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------------------