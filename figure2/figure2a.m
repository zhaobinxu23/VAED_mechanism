%% nature_Extended_data_figure1_A_beautified.m
% Beautified version for Extended Data Figure 1A
% Left semicircle: k_on
% Right semicircle: k_off

clear; clc; close all;

%% ==========================================================
% 1. Figure setup
% ==========================================================
fig = figure('Color', 'w', ...
             'Position', [100, 100, 1500, 850], ...
             'Renderer', 'painters');

ax_main = axes('Position', [0.06, 0.13, 0.74, 0.76]);
hold(ax_main, 'on');
axis(ax_main, 'equal');

set(ax_main, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 1.0, ...
    'TickDir', 'out', ...
    'Layer', 'top', ...
    'Color', 'none');

%% ==========================================================
% 2. Parameters and colormaps
% ==========================================================

% k_on range: 5 x 10^-22 to 5 x 10^-13
min_log_kon = log10(5 * 10^(-22));
max_log_kon = log10(5 * 10^(-13));

% k_off range: 10^0 to 10^9
min_log_koff = 0;
max_log_koff = 9;

% Nature-style soft blue colormap
blues_map = interp1([1 128 256], ...
    [0.92 0.96 1.00; ...   % very light blue
     0.35 0.65 0.90; ...   % medium blue
     0.05 0.20 0.55], ...  % deep blue
    1:256);

% Nature-style soft red colormap
reds_map = interp1([1 128 256], ...
    [1.00 0.94 0.92; ...   % very light red
     0.95 0.45 0.35; ...   % medium red
     0.55 0.05 0.08], ...  % deep red
    1:256);

% Half-circle angles
theta_left  = linspace(pi/2, 3*pi/2, 80);
theta_right = linspace(-pi/2, pi/2, 80);

%% ==========================================================
% 3. Generate data
% ==========================================================

data_struct = [];

for m = 1:10
    for n = 1:10
        ab_id = 10 * (m - 1) + n;
        
        val_kon = 5 * 10^(m - 23);
        val_koff = 10^(n - 1);
        
        log_kon_val = log10(val_kon);
        log_koff_val = log10(val_koff);
        
        % Affinity group coordinate
        x_val = -(n - m) - 22;
        
        data_struct = [data_struct; ...
            ab_id, x_val, log_kon_val, log_koff_val, m]; %#ok<AGROW>
    end
end

% Sort by affinity group and m
[~, sort_idx] = sortrows(data_struct, [2, 5]);
data_struct = data_struct(sort_idx, :);

%% ==========================================================
% 4. Calculate node positions
% ==========================================================

radius = 0.35;
x_counts = containers.Map('KeyType', 'double', 'ValueType', 'double');
x_totals = containers.Map('KeyType', 'double', 'ValueType', 'double');

% Count number of antibodies in each x group
for i = 1:size(data_struct, 1)
    x_val = data_struct(i, 2);
    
    if isKey(x_totals, x_val)
        x_totals(x_val) = x_totals(x_val) + 1;
    else
        x_totals(x_val) = 1;
    end
end

unique_x = sort(cell2mat(keys(x_totals)));

%% ==========================================================
% 5. Light guide lines behind nodes
% ==========================================================

xlim(ax_main, [-32.8, -11.2]);
ylim(ax_main, [-6.4, 6.2]);

for xx = unique_x
    plot(ax_main, [xx xx], ylim(ax_main), ...
        '-', 'Color', [0.92 0.92 0.92], 'LineWidth', 0.6);
end

%% ==========================================================
% 6. Draw nodes
% ==========================================================

for i = 1:size(data_struct, 1)
    id     = data_struct(i, 1);
    x      = data_struct(i, 2);
    l_kon  = data_struct(i, 3);
    l_koff = data_struct(i, 4);
    
    % Count current position in this x group
    if isKey(x_counts, x)
        x_counts(x) = x_counts(x) + 1;
    else
        x_counts(x) = 1;
    end
    
    count_now = x_counts(x);
    total_now = x_totals(x);
    
    % Vertically centered stacking
    y_spacing = 0.82;
    y = (count_now - (total_now + 1) / 2) * y_spacing;
    
    % ---------- Left semicircle: k_on ----------
    norm_idx = floor((l_kon - min_log_kon) / ...
        (max_log_kon - min_log_kon) * 255) + 1;
    norm_idx = max(1, min(256, norm_idx));
    c_left = blues_map(norm_idx, :);
    
    X_L = x + radius * cos(theta_left);
    Y_L = y + radius * sin(theta_left);
    
    patch(ax_main, ...
        'XData', X_L, ...
        'YData', Y_L, ...
        'FaceColor', c_left, ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.45);
    
    % ---------- Right semicircle: k_off ----------
    norm_idx_off = floor((l_koff - min_log_koff) / ...
        (max_log_koff - min_log_koff) * 255) + 1;
    norm_idx_off = max(1, min(256, norm_idx_off));
    c_right = reds_map(norm_idx_off, :);
    
    X_R = x + radius * cos(theta_right);
    Y_R = y + radius * sin(theta_right);
    
    patch(ax_main, ...
        'XData', X_R, ...
        'YData', Y_R, ...
        'FaceColor', c_right, ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.45);
    
    % ---------- Text label ----------
    if norm_idx > 155 || norm_idx_off > 155
        txt_col = [1 1 1];
    else
        txt_col = [0.12 0.12 0.12];
    end
    
    text(ax_main, x, y, sprintf('%d', id), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 7.5, ...
        'FontWeight', 'bold', ...
        'Color', txt_col);
end

%% ==========================================================
% 7. Main axis decoration
% ==========================================================

xlim(ax_main, [-32.8, -11.2]);
ylim(ax_main, [-6.4, 6.2]);
xticks(ax_main, unique_x);

xlabel(ax_main, 'Affinity group, log_{10}(K_D)', ...
    'FontName', 'Arial', ...
    'FontSize', 13, ...
    'FontWeight', 'bold');

set(ax_main, ...
    'YTick', [], ...
    'YColor', 'none', ...
    'Box', 'off', ...
    'XColor', [0.15 0.15 0.15], ...
    'FontName', 'Arial', ...
    'FontSize', 11);

text(ax_main, -31.0, -6.75, '\leftarrow weaker binding', ...
    'Color', [0.55 0.10 0.10], ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

text(ax_main, -13.0, -6.75, 'stronger binding \rightarrow', ...
    'Color', [0.05 0.38 0.15], ...
    'FontName', 'Arial', ...
    'FontSize', 12, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

title(ax_main, 'IgM affinity clusters mapped by kinetic rates', ...
    'FontName', 'Arial', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

%% ==========================================================
% 8. Colorbars
% ==========================================================

% ---------- k_on colorbar ----------
ax_cb1 = axes('Position', [0.84, 0.18, 0.025, 0.30], 'Visible', 'off');
colormap(ax_cb1, blues_map);
caxis(ax_cb1, [min_log_kon max_log_kon]);

cb1 = colorbar(ax_cb1, 'Position', [0.84, 0.18, 0.025, 0.30]);
cb1.Label.String = 'k_{on} value';
cb1.Label.FontName = 'Arial';
cb1.Label.FontSize = 11;
cb1.Label.FontWeight = 'bold';
cb1.FontName = 'Arial';
cb1.FontSize = 9;
cb1.LineWidth = 0.8;
cb1.Box = 'off';

% Ticks for k_on
kon_exponents = -22:-13;
kon_tick_values = log10(5 * 10.^kon_exponents);
cb1.Ticks = kon_tick_values;

kon_labels = cell(numel(kon_exponents), 1);
for k = 1:numel(kon_exponents)
    kon_labels{k} = sprintf('5×10^{%d}', kon_exponents(k));
end
cb1.TickLabels = kon_labels;

% ---------- k_off colorbar ----------
ax_cb2 = axes('Position', [0.91, 0.18, 0.025, 0.30], 'Visible', 'off');
colormap(ax_cb2, reds_map);
caxis(ax_cb2, [min_log_koff max_log_koff]);

cb2 = colorbar(ax_cb2, 'Position', [0.91, 0.18, 0.025, 0.30]);
cb2.Label.String = 'k_{off} value';
cb2.Label.FontName = 'Arial';
cb2.Label.FontSize = 11;
cb2.Label.FontWeight = 'bold';
cb2.FontName = 'Arial';
cb2.FontSize = 9;
cb2.LineWidth = 0.8;
cb2.Box = 'off';

% Ticks for k_off
koff_exponents = 0:9;
cb2.Ticks = koff_exponents;

koff_labels = cell(numel(koff_exponents), 1);
for k = 1:numel(koff_exponents)
    koff_labels{k} = sprintf('10^{%d}', koff_exponents(k));
end
cb2.TickLabels = koff_labels;

%% ==========================================================
% 9. Legend inset
% ==========================================================

ax_legend = axes('Position', [0.805, 0.72, 0.18, 0.18], 'Visible', 'off');
hold(ax_legend, 'on');
xlim(ax_legend, [0 1]);
ylim(ax_legend, [0 1]);
axis(ax_legend, 'equal');

% Background panel
rectangle(ax_legend, ...
    'Position', [0.02, 0.05, 0.96, 0.90], ...
    'Curvature', 0.08, ...
    'FaceColor', [0.98 0.98 0.98], ...
    'EdgeColor', [0.82 0.82 0.82], ...
    'LineWidth', 0.8);

% Demo node
demo_r = 0.16;
dx = 0.33;
dy = 0.56;

patch('Parent', ax_legend, ...
      'XData', dx + demo_r * cos(theta_left), ...
      'YData', dy + demo_r * sin(theta_left), ...
      'FaceColor', blues_map(190, :), ...
      'EdgeColor', [1 1 1], ...
      'LineWidth', 0.5);

patch('Parent', ax_legend, ...
      'XData', dx + demo_r * cos(theta_right), ...
      'YData', dy + demo_r * sin(theta_right), ...
      'FaceColor', reds_map(190, :), ...
      'EdgeColor', [1 1 1], ...
      'LineWidth', 0.5);

text(ax_legend, 0.58, 0.66, 'Node encoding', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 10, ...
    'HorizontalAlignment', 'left');

text(ax_legend, 0.50, 0.34, 'left: k_{on}    right: k_{off}', ...
    'FontName', 'Arial', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 8.5, ...
    'Color', [0.25 0.25 0.25]);

%% ==========================================================
% 10. Export high-resolution figure
% ==========================================================

set(fig, 'PaperPositionMode', 'auto');

exportgraphics(fig, 'nature_Extended_data_figure1_A_beautified.png', ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

exportgraphics(fig, 'nature_Extended_data_figure1_A_beautified.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

fprintf('Figure exported successfully:\n');
fprintf('  - nature_Extended_data_figure1_A_beautified.png\n');
fprintf('  - nature_Extended_data_figure1_A_beautified.pdf\n');