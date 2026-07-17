%% nature_Extended_data_figure1_B_beautified.m
% Beautified Extended Data Figure 1B
% AA distribution generated from normal CDF with mu = 9, sigma = 0.8

clear; clc; close all;

%% ==========================================================
% 1. Figure setup
% ==========================================================
fig = figure('Color', 'w', ...
             'Position', [100, 100, 1500, 850], ...
             'Renderer', 'painters');

ax_main = axes('Position', [0.07, 0.14, 0.72, 0.76]);
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
% 2. Data preparation
% ==========================================================

% ----------------------------------------------------------
% Generate AA_raw distribution from normal CDF
% This replaces the original hard-coded raw_str values
% ----------------------------------------------------------
mu = 9;        % 均值
sigma = 0.8;   % 标准差

prob_A = zeros(1, 10);

% prob_A(1) = 0.001;
% prob_A(2:9) = 0;
prob_A(1) = normcdf(5, mu, sigma);

% 计算累积分布概率
for i = 2:5
    prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
end

% 对称分布
for i = 6:10
    prob_A(i) = prob_A(11-i);
end

AA_raw = zeros(1, 100);
for i = 1:10
    for j = 1:10
        AA_raw(10*(i-1)+j) = prob_A(i) * prob_A(j);
    end
end

% 如果你仍然需要 raw_str，这里自动生成一个字符串版本
raw_str = sprintf('%.16g ', AA_raw);

% 绝对数量计算
TOTAL_POOL = 4 * 10^19;
AA = AA_raw * TOTAL_POOL;

%% ==========================================================
% 3. Build plotting data
% ==========================================================

plot_data = [];
for m = 1:10
    for n = 1:10
        idx = 10*(m-1) + n;
        val = AA(idx);
        x_val = -(n - m) - 22;
        plot_data = [plot_data; idx, x_val, val, m, n]; %#ok<AGROW>
    end
end

% 排序保证每个 affinity group 内部堆叠稳定
[~, sort_idx] = sortrows(plot_data, [2, 4, 5]);
plot_data = plot_data(sort_idx, :);

unique_x = unique(plot_data(:, 2));
unique_x = sort(unique_x);

%% ==========================================================
% 4. Colormap and value scaling
% ==========================================================

% Nature-like purple-orange colormap
% Low values: light grey-purple
% High values: deep orange-red
cmap = interp1([1 80 160 256], ...
    [0.94 0.94 0.97; ...   % very light grey-purple
     0.68 0.76 0.88; ...   % soft blue
     0.98 0.70 0.42; ...   % warm orange
     0.70 0.10 0.10], ...  % deep red
    1:256);

colormap(ax_main, cmap);

% 避免 log10(0)
positive_vals = plot_data(plot_data(:,3) > 0, 3);
min_val = min(positive_vals);
max_val = max(positive_vals);

min_log = log10(min_val);
max_log = log10(max_val);

%% ==========================================================
% 5. Light guide lines
% ==========================================================

xlim(ax_main, [-32.8, -11.2]);
ylim(ax_main, [-6.4, 6.2]);

for xx = unique_x'
    plot(ax_main, [xx xx], ylim(ax_main), ...
        '-', 'Color', [0.92 0.92 0.92], 'LineWidth', 0.6);
end

%% ==========================================================
% 6. Draw nodes
% ==========================================================

radius = 0.35;
theta = linspace(0, 2*pi, 120);

count_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
total_map = containers.Map('KeyType', 'double', 'ValueType', 'double');

% 统计每列总数
for i = 1:size(plot_data, 1)
    x = plot_data(i, 2);
    if isKey(total_map, x)
        total_map(x) = total_map(x) + 1;
    else
        total_map(x) = 1;
        count_map(x) = 0;
    end
end

for i = 1:size(plot_data, 1)
    id  = plot_data(i, 1);
    x   = plot_data(i, 2);
    val = plot_data(i, 3);
    m   = plot_data(i, 4);
    n   = plot_data(i, 5);

    % Y 坐标：居中堆叠
    count_map(x) = count_map(x) + 1;
    y_spacing = 0.82;
    y = (count_map(x) - (total_map(x) + 1) / 2) * y_spacing;

    % 颜色映射，使用 log10 数值
    if val <= 0
        color_idx = 1;
    else
        val_log = log10(val);
        color_idx = floor((val_log - min_log) / (max_log - min_log) * 255) + 1;
        color_idx = max(1, min(256, color_idx));
    end

    face_col = cmap(color_idx, :);

    % 画圆点
    X = x + radius * cos(theta);
    Y = y + radius * sin(theta);

    patch(ax_main, ...
        'XData', X, ...
        'YData', Y, ...
        'FaceColor', face_col, ...
        'EdgeColor', [1 1 1], ...
        'LineWidth', 0.55);

    % 根据颜色深浅自动调整文字颜色
    brightness = 0.299 * face_col(1) + 0.587 * face_col(2) + 0.114 * face_col(3);
    if brightness < 0.55
        txt_col = [1 1 1];
    else
        txt_col = [0.12 0.12 0.12];
    end

    % 标签：默认显示编号
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

title(ax_main, 'IgM repertoire abundance landscape', ...
    'FontName', 'Arial', ...
    'FontSize', 16, ...
    'FontWeight', 'bold');

%% ==========================================================
% 8. Colorbar
% ==========================================================

ax_cb = axes('Position', [0.84, 0.18, 0.035, 0.50], 'Visible', 'off');
colormap(ax_cb, cmap);
caxis(ax_cb, [min_log max_log]);

cb = colorbar(ax_cb, 'Position', [0.84, 0.18, 0.035, 0.50]);
cb.Label.String = 'Absolute abundance';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';
cb.FontName = 'Arial';
cb.FontSize = 9;
cb.LineWidth = 0.8;
cb.Box = 'off';

% 设置 colorbar ticks 为 10 的幂次
log_ticks = ceil(min_log):floor(max_log);
cb.Ticks = log_ticks;

cb_labels = cell(numel(log_ticks), 1);
for k = 1:numel(log_ticks)
    cb_labels{k} = sprintf('10^{%d}', log_ticks(k));
end
cb.TickLabels = cb_labels;

%% ==========================================================
% 9. Distribution inset: show prob_A
% ==========================================================

ax_inset = axes('Position', [0.805, 0.73, 0.18, 0.18]);
hold(ax_inset, 'on');

bar(ax_inset, 1:10, prob_A, ...
    'FaceColor', [0.35 0.55 0.78], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.85);

plot(ax_inset, 1:10, prob_A, '-o', ...
    'Color', [0.12 0.25 0.45], ...
    'MarkerFaceColor', [1 1 1], ...
    'MarkerSize', 4, ...
    'LineWidth', 1.2);

set(ax_inset, ...
    'FontName', 'Arial', ...
    'FontSize', 8.5, ...
    'LineWidth', 0.8, ...
    'TickDir', 'out', ...
    'Box', 'off', ...
    'XColor', [0.2 0.2 0.2], ...
    'YColor', [0.2 0.2 0.2]);

xlabel(ax_inset, 'Index', 'FontName', 'Arial', 'FontSize', 8.5);
ylabel(ax_inset, 'P', 'FontName', 'Arial', 'FontSize', 8.5);
title(ax_inset, sprintf('Normal CDF\n\mu = %.1f, \sigma = %.1f', mu, sigma), ...
    'FontName', 'Arial', ...
    'FontSize', 9, ...
    'FontWeight', 'bold');

xlim(ax_inset, [0.5 10.5]);

%% ==========================================================
% 10. Optional annotation
% ==========================================================

annotation(fig, 'textbox', [0.805, 0.065, 0.18, 0.08], ...
    'String', sprintf('Total pool = %.1e\nAA_{ij} = P_i P_j \times total pool', TOTAL_POOL), ...
    'FontName', 'Arial', ...
    'FontSize', 9, ...
    'Color', [0.25 0.25 0.25], ...
    'EdgeColor', [0.85 0.85 0.85], ...
    'BackgroundColor', [0.98 0.98 0.98], ...
    'FitBoxToText', 'on');

%% ==========================================================
% 11. Export high-resolution figure
% ==========================================================

set(fig, 'PaperPositionMode', 'auto');

exportgraphics(fig, 'nature_Extended_data_figure1_B_beautified.png', ...
    'Resolution', 600, ...
    'BackgroundColor', 'white');

exportgraphics(fig, 'nature_Extended_data_figure1_B_beautified.pdf', ...
    'ContentType', 'vector', ...
    'BackgroundColor', 'white');

%% ==========================================================
% 12. Print generated distribution
% ==========================================================

fprintf('\nGenerated prob_A distribution:\n');
disp(prob_A);

fprintf('Generated AA_raw distribution:\n');
disp(AA_raw);

fprintf('Generated raw_str:\n');
fprintf('%s\n', raw_str);

fprintf('\nFigure exported successfully:\n');
fprintf('  - nature_Extended_data_figure1_B_beautified.png\n');
fprintf('  - nature_Extended_data_figure1_B_beautified.pdf\n');