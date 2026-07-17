%  0. Uniform time grid and basic indices
%  ------------------------------------------------------------
load('E:\antibody_dynamics\epitope_competition\figure6\batch_results_vaccine_variant\challenge_alpha_0.mat') %% change it into the corresponding location of challenge_alpha_0-4.mat
t_grid = (0:1:1000)';

N = p_chal.N;
age_vec = p_chal.a_vec(:);

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
Y_first_grid = interp1(t_chal, y_chal, t_grid, 'linear', 'extrap');
% Y_second_grid = interp1(t1, y1, t_grid, 'linear', 'extrap');

Y_first_grid = real(Y_first_grid);
% Y_second_grid = real(Y_second_grid);

Y_first_grid(~isfinite(Y_first_grid)) = 0;
% Y_second_grid(~isfinite(Y_second_grid)) = 0;

% 数值误差可能产生极小负数，统一截断
Y_first_grid(Y_first_grid < 0) = 0;
% Y_second_grid(Y_second_grid < 0) = 0;

%% ------------------------------------------------------------
%  2. Extract variables on t_grid
%  ------------------------------------------------------------

V_first  = Y_first_grid(:, idx_V);
% V_second = Y_second_grid(:, idx_V);

% infectious virus = free virus + X(3603)
infectious_V_first  = Y_first_grid(:, idx_V)  + Y_first_grid(:, idx_X3603);
% infectious_V_second = Y_second_grid(:, idx_V) + Y_second_grid(:, idx_X3603);

I_first  = Y_first_grid(:, idx_I);
% I_second = Y_second_grid(:, idx_I);

total_I_first  = sum(I_first, 2);
% total_I_second = sum(I_second, 2);

cum_nat_first  = Y_first_grid(:, idx_cum_nat);
cum_adcc_first = Y_first_grid(:, idx_cum_adcc);
cum_tc_first   = Y_first_grid(:, idx_cum_tc);

% cum_nat_second  = Y_second_grid(:, idx_cum_nat);
% cum_adcc_second = Y_second_grid(:, idx_cum_adcc);
% cum_tc_second   = Y_second_grid(:, idx_cum_tc);

cum_total_first  = cum_nat_first  + cum_adcc_first  + cum_tc_first;
% cum_total_second = cum_nat_second + cum_adcc_second + cum_tc_second;

%% ============================================================
%  Figure 4a
%  Free extracellular virus and infectious virus
%  ============================================================

figure('Color','w','Position',[100 100 1150 480]);

subplot(1,2,1);
hold on;
plot(t_grid, V_first, 'k-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Free extracellular virus');
legend({'mutant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

subplot(1,2,2);
hold on;
plot(t_grid, infectious_V_first, 'k-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Infectious virus');

legend({'mutant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

sgtitle('Extracellular viral load during secondary infection');

saveas(gcf, 'figure4a_viral_load_free_and_infectious.png');
saveas(gcf, 'figure4a_viral_load_free_and_infectious.fig');

%% ============================================================
%  Figure 4b
%  Total infected-cell burden over time
%  ============================================================

figure('Color','w','Position',[100 100 650 500]);
hold on;
plot(t_grid, total_I_first, 'k-', 'LineWidth', 2);
set(gca, 'YScale', 'log');
xlabel('Time after infection');
ylabel('Total infected cells');
title('Total infected-cell burden over time');
legend({'mutant infection'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);

saveas(gcf, 'figure4b_total_infected_cell_burden.png');
saveas(gcf, 'figure4b_total_infected_cell_burden.fig');

%% ============================================================
%  Figure 4c
%  Cumulative infected-cell lysis
%  ============================================================

figure('Color','w','Position',[100 100 1200 520]);


hold on;
plot(t_grid, cum_nat_first,   'b-',  'LineWidth', 2);
plot(t_grid, cum_adcc_first,  'r-',  'LineWidth', 2);
plot(t_grid, cum_tc_first,    'g-',  'LineWidth', 2);
plot(t_grid, cum_total_first, 'k--', 'LineWidth', 2);
xlabel('Time after infection');
ylabel('Cumulative lysed infected cells');
legend({'Natural lysis','ADCC lysis','Tc killing','Total'}, 'Location','best');
box on;
set(gca, 'FontSize', 12);


sgtitle('Cumulative infected-cell lysis');

saveas(gcf, 'figure4c_cumulative_infected_cell_lysis.png');
saveas(gcf, 'figure4c_cumulative_infected_cell_lysis.fig');


%% ============================================================
%  Figure 4d
%  Distribution of infected-cell by infection age
%  ============================================================

Z_I_first  = log10(I_first' + 1);


Z_I_first(~isfinite(Z_I_first)) = 0;


common_I_min = 0;
common_I_max = max(Z_I_first(:));
if ~isfinite(common_I_max) || common_I_max <= common_I_min
    common_I_max = 1;
end

figure('Color','w','Position',[100 100 1200 500]);

imagesc(t_grid, age_vec, Z_I_first);
set(gca, 'YDir', 'normal');
caxis([common_I_min common_I_max]);
colormap hot;
cb1 = colorbar;
ylabel(cb1, 'log_{10}(infected cells + 1)');
xlabel('Time after infection');
ylabel('Infection age');
set(gca, 'FontSize', 12);



sgtitle('Distribution of infected cells by infection age');

saveas(gcf, 'figure4d_infected_cell_distribution_by_age.png');
saveas(gcf, 'figure4d_infected_cell_distribution_by_age.fig');


%% ============================================================
%  Combined Figure 4
%  Viral load, infected-cell burden, cell lysis and
%  infected-cell age distribution
%  ============================================================

% -------------------------------------------------------------
% Prepare heat-map data
% --------------------------------------------------------------
Z_I_first = log10(I_first' + 1);
Z_I_first(~isfinite(Z_I_first)) = 0;

common_I_min = 0;
common_I_max = max(Z_I_first(:));

if ~isfinite(common_I_max) || common_I_max <= common_I_min
    common_I_max = 1;
end

% -------------------------------------------------------------
% Create a single combined figure
% -------------------------------------------------------------
fig4 = figure( ...
    'Color', 'w', ...
    'Position', [80 50 1300 1050]);

% Three rows and two columns:
% Row 1: extracellular virus and infectious virus
% Row 2: infected-cell burden and cumulative cell lysis
% Row 3: infected-cell age distribution spanning two columns
tl = tiledlayout(fig4, 3, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% ------------------------------------------------------------
% Free extracellular virus
% -------------------------------------------------------------
ax1 = nexttile(tl, 1);

plot(ax1, t_grid, V_first, ...
    'k-', ...
    'LineWidth', 2);

set(ax1, ...
    'YScale', 'log', ...
    'FontSize', 12, ...
    'LineWidth', 1, ...
    'Box', 'on');

xlabel(ax1, 'Time after infection');
ylabel(ax1, 'Free extracellular virus');
title(ax1, 'Extracellular viral load');

legend(ax1, {'Mutant infection'}, ...
    'Location', 'best', ...
    'Box', 'off');

%% ------------------------------------------------------------
% Infectious virus
% -------------------------------------------------------------
ax2 = nexttile(tl, 2);

plot(ax2, t_grid, infectious_V_first, ...
    'k-', ...
    'LineWidth', 2);

set(ax2, ...
    'YScale', 'log', ...
    'FontSize', 12, ...
    'LineWidth', 1, ...
    'Box', 'on');

xlabel(ax2, 'Time after infection');
ylabel(ax2, 'Infectious virus');
title(ax2, 'Infectious viral load');

legend(ax2, {'Mutant infection'}, ...
    'Location', 'best', ...
    'Box', 'off');

%% ------------------------------------------------------------
% Total infected-cell burden
% -------------------------------------------------------------
ax3 = nexttile(tl, 3);

plot(ax3, t_grid, total_I_first, ...
    'k-', ...
    'LineWidth', 2);

set(ax3, ...
    'YScale', 'log', ...
    'FontSize', 12, ...
    'LineWidth', 1, ...
    'Box', 'on');

xlabel(ax3, 'Time after infection');
ylabel(ax3, 'Total infected cells');
title(ax3, 'Total infected-cell burden');

legend(ax3, {'Mutant infection'}, ...
    'Location', 'best', ...
    'Box', 'off');

%% ------------------------------------------------------------
% Cumulative infected-cell lysis
% -------------------------------------------------------------
ax4 = nexttile(tl, 4);

hold(ax4, 'on');

plot(ax4, t_grid, cum_nat_first, ...
    'b-', ...
    'LineWidth', 2);

plot(ax4, t_grid, cum_adcc_first, ...
    'r-', ...
    'LineWidth', 2);

plot(ax4, t_grid, cum_tc_first, ...
    'Color', [0.10 0.55 0.10], ...
    'LineStyle', '-', ...
    'LineWidth', 2);

plot(ax4, t_grid, cum_total_first, ...
    'k--', ...
    'LineWidth', 2);

hold(ax4, 'off');

set(ax4, ...
    'FontSize', 12, ...
    'LineWidth', 1, ...
    'Box', 'on');

xlabel(ax4, 'Time after infection');
ylabel(ax4, 'Cumulative lysed infected cells');
title(ax4, 'Cumulative infected-cell lysis');

legend(ax4, ...
    {'Natural lysis', 'ADCC lysis', 'Tc killing', 'Total'}, ...
    'Location', 'best', ...
    'Box', 'off');

%% ------------------------------------------------------------
% Distribution of infected cells by infection age
% This panel spans the two columns of the final row
% -------------------------------------------------------------
ax5 = nexttile(tl, 5, [1 2]);

imagesc(ax5, t_grid, age_vec, Z_I_first);

set(ax5, ...
    'YDir', 'normal', ...
    'FontSize', 12, ...
    'LineWidth', 1, ...
    'Box', 'on');

caxis(ax5, [common_I_min, common_I_max]);
colormap(ax5, hot);

cb1 = colorbar(ax5);
ylabel(cb1, 'log_{10}(infected cells + 1)', ...
    'FontSize', 12);

xlabel(ax5, 'Time after infection');
ylabel(ax5, 'Infection age');
title(ax5, 'Distribution of infected cells by infection age');

%% ------------------------------------------------------------
% Overall title
% Delete this line if the journal figure should have no title
% -------------------------------------------------------------
sgtitle(tl, ...
    'Dynamics of viral infection and infected-cell clearance', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

%% ------------------------------------------------------------
% Save the combined figure
% -------------------------------------------------------------
savefig(fig4, 'figure4_combined.fig');

% Export high-resolution raster image
exportgraphics(fig4, ...
    'figure4_combined.png', ...
    'Resolution', 600);

% Export vector PDF, suitable for manuscript preparation
exportgraphics(fig4, ...
    'figure4_combined.pdf', ...
    'ContentType', 'vector');

  