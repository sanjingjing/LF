clear all;
clc; 
close all

% ===================== 1. 定义模型常数 =====================
t1 = 600;      % 温度增量 (K)
t2 = 1200;     % 环境温度 (K)
m = 3;       
n = 0.9;     
R = 0.04;
Z = 0.4;
x = 0;

[y,z] = meshgrid(linspace(-0.04, 0.04, 500), ...
                 linspace(0, 0.4, 500));

term = m*(x.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;

% 原始物理温度场
T = t1 * exp(-term.^2) + t2;

% ===================== 2. 只用于绘图的温度场 =====================
T_plot = T;

% 注意：因为环境温度 t2 = 1200 K
% 所以 T_cut 必须大于 1200，火焰外区域才会被隐藏
T_cut = 1250;    

T_plot(T <= T_cut) = NaN;

% ===================== 3. 绘图 =====================
figure('Color','white');
ax = axes;
ax.Color = 'white';

h = pcolor(y,z,T_plot);
shading interp;
set(h,'EdgeColor','none');

colormap(jet(256));

% 固定色标范围，与论文一致
clim([1200 1800]);

cb = colorbar;
cb.Ticks = 1200:75:1800;
cb.Label.String = '温度 (K)';
cb.Label.Rotation = 0;
cb.Label.FontSize = 12;

xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');

pbaspect([1 5 1]);
box on;
grid off;