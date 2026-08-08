clear all;
clc;
close all;

%% ===================== 1. 模型参数 =====================
t1 = 1200;      % 温度增量 (K)
t2 = 900;       % 环境温度下限 (K)
m  = 10;
n  = 0.4;
R  = 0.025;     % 火焰半径 (m)
Z  = 0.075;     % 火焰高度尺度 (m)

% 固定截面位置
x = 0;          % m

% 为了显示风格接近论文，横向和纵向范围稍微放大一点
y_min = -0.05;  % m
y_max =  0.05;  % m
z_min = -0.015; % m
z_max =  0.085; % m

Ny = 1000;
Nz = 1000;

[y, z] = meshgrid(linspace(y_min, y_max, Ny), ...
                  linspace(z_min, z_max, Nz));

%% ===================== 2. 计算温度场 =====================
T = NaN(size(y));

% 只在 z>=0 的物理区域内计算
mask_phys = (z >= 0);

term = m .* ( z.^2 ./ Z.^2 + (x.^2 + y.^2) ./ R.^2 ) - n;
T_all = t1 .* exp(-(term.^2)) + t2;

T(mask_phys) = T_all(mask_phys);

%% ===================== 3. 仅用于显示的温度场 =====================
T_plot = T;

% 截断低温区，让火焰边界更像论文图
T_cut = 920;
T_plot(T_plot <= T_cut) = NaN;

%% ===================== 4. 自定义颜色设置 =====================
% 深红（低温） -> 红 -> 橙 -> 黄（高温）
base_colors = [
    0.25 0.00 0.00   % 深红：低温
    0.60 0.00 0.00   % 红
    0.90 0.00 0.00   % 亮红
    1.00 0.35 0.00   % 橙红
    1.00 0.65 0.00   % 橙
    1.00 0.85 0.00   % 金黄
    1.00 1.00 0.00   % 黄：高温
];

x_base = linspace(0,1,size(base_colors,1));
xq = linspace(0,1,256);
cmap = interp1(x_base, base_colors, xq);

% 温度显示范围
Tmin_show = 900;
Tmax_show = 2300;

%% ===================== 5. 绘图 =====================
figure('Color', 'white', 'Position', [200 100 520 420]);

ax = axes;
hold(ax, 'on');
ax.Color = 'black';   % NaN区域显示为黑色

h = pcolor(ax, y*1000, z*1000, T_plot);   % 转成 mm
shading(ax, 'interp');
set(h, 'EdgeColor', 'none');

colormap(ax, cmap);
clim(ax, [Tmin_show Tmax_show]);

%% ===================== 6. 坐标轴样式 =====================
xlim(ax, [-50 50]);
ylim(ax, [-15 85]);

xlabel(ax, 'Width [mm]', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
ylabel(ax, 'Height [mm]', ...
    'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');

set(ax, ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 1.2, ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'Box', 'on');

xticks(ax, [-50 -25 0 25 50]);
yticks(ax, [0 25 50 75]);

% 让图像比例接近论文
pbaspect(ax, [1 1 1.05]);

%% ===================== 7. 图内文字 =====================
text(ax, -47, 84, '(a)', ...
    'Color', 'w', ...
    'FontSize', 28, ...
    'FontWeight', 'bold', ...
    'VerticalAlignment', 'top');

text(ax, 0, 84, 'Original', ...
    'Color', 'w', ...
    'FontSize', 20, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top');

%% ===================== 8. colorbar =====================
cb = colorbar(ax);
cb.Position = [0.84 0.16 0.05 0.67];
cb.FontSize = 12;
cb.FontWeight = 'bold';
cb.Ticks = [1000 1200 1400 1600 1800 2000 2200];
cb.Label.String = 'Temperature [K]';
cb.Label.FontSize = 16;
cb.Label.FontWeight = 'bold';

%% ===================== 9. 保存图片 =====================
exportgraphics(gcf, 'flame_original_style_red_to_yellow.png', 'Resolution', 300);