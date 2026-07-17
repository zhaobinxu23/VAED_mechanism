function [p, y0, t_span] = build_primary_vaccine_model(vaccine_type, dose)

%% =========================
%  1. 参数设置
% =========================

p.N  = 400;
p.da = 1.0;
p.dt_step = 1;
p.infection_off = false;
p.extinction_threshold = 1;

%% 免疫系统参数
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
%% =========================
%  2. 疫苗接种阶段：禁止病毒入侵
% =========================

p.k4 = 0;          % 关键：疫苗接种阶段病毒不入侵细胞
p.km = 1.0e6;
p.k6 = 1.0e8;
p.k7 = 0.01;
p.c_clear = 0.1;
p.Tc_generation = 5e-2;

p.v_start = 0;     % 疫苗阶段不产生胞内病毒
p.k5 = 0;
p.a_vec = (0:p.N-1)' * p.da;
p.vin_vec = zeros(p.N, 1);   % 关键：不同感染龄胞内病毒浓度均为0

p.Tc_binding = 1e-5;
p.theta_lysis = 1e20;
p.theta_adcc  = 1e10;
p.theta_tc = 1e8;
p.n_hill = 2;
p.k_lysis_max = 1;
p.k_adcc_max = 1;
p.k_tc_max = 1;
p.k_kill_tc = 0.1;

%% =========================
%  3. 标记疫苗类型
% =========================

p.vaccine_type = vaccine_type;

switch vaccine_type
    case 'nonRBD'
        % non-RBD 区域接种：
        % RBD 抗体 1801:3600 不参与病毒结合
        p.binding_mode = 'nonRBD_only';

    case 'RBD'
        % RBD 区域接种：
        % non-RBD 抗体 1:1800 不参与病毒结合
        p.binding_mode = 'RBD_only';

    case 'inactivated'
        % 全病毒灭活疫苗：
        % non-RBD 和 RBD 都参与结合
        p.binding_mode = 'both';

    otherwise
        error('Unknown vaccine_type: %s', vaccine_type);
end

%% =========================
%  4. 初始条件
% =========================

T0 = 1e10;
I0 = zeros(p.N, 1);  % 关键：疫苗阶段无感染细胞

% 疫苗接种阶段，原则上不应该有真正可复制的感染病毒。
% 如果你的模型中 V 表示可入侵病毒，建议设为0。
V0 = dose;

Tc0 = 1e5;

x0 = zeros(3604, 1);

%% ----------------------------------------------------------
%  抗体初始分布部分
%  这里建议直接复制你原始代码中关于 x0(1:3600) 初始化的部分
% ----------------------------------------------------------
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
% 示例：
% mu = 9;
% sigma = 0.80;
% ...
% x0(1:3600) = ...

%% ----------------------------------------------------------
%  疫苗抗原剂量
% ----------------------------------------------------------
% 根据你原始代码：
% x0(3602) = 1e1; % vaccine or virus
%
% 因此建议把疫苗剂量放在 x0(3602)。
% 如果你的模型实际上使用 V0 代表疫苗抗原，则需要改成 V0 = dose。
% 但从你附件片段看，x0(3602) 更像是 vaccine or virus 变量。

x0(3601) = 0;       % 如果原始代码中这里是 E，请替换为你的 E
x0(3602) = dose;    % 疫苗剂量
x0(3603) = 0;       % V-A1 complex
x0(3604) = 0;       % V-A2 complex

p.AA = AA;
p.AA_new = AA_new;


cum_death_0 = [0; 0; 0];
L_age_0 = zeros(p.N, 1);

y0 = [T0; ...
      I0; ...
      V0; ...
      Tc0; ...
      x0; ...
      cum_death_0; ...
      L_age_0];

%% =========================
%  5. 模拟时间
% =========================

t_span = [0, 1000];

end