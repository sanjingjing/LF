clear;
clc;

% ===================== 1. 定义模型常数 =====================
t1 = 1200;     % 温度增量 (K)
t2 = 900;      % 环境温度 (K)
m  = 3;      
n  = 0.9;      % 形状参数
R  = 0.04;
Z  = 0.4;
x  = 0;     % 固定切面位置

% 网格
Ny = 500;
Nz = 500;

[y,z] = meshgrid(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz));

% ===================== 2. 温度场 =====================
term = m * (x.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
T = t1 * exp(-(term.^2)) + t2;

% ===================== 3. 单色辐射模型（论文式 2-13） =====================
lambda = 0.610;          % 波长 (um)
c1 = 3.7418e-16;         % 第一辐射常数
c2 = 1.4388e4;           % 第二辐射常数 (um*K)

% 单色黑体辐射
Eb = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T)) - 1);

% ===================== 4. 辐射传输衰减（Beer 定律） =====================
% 这里把当前二维图理解为：沿 y 方向传播到相机
% 每个像素到相机的传播距离近似取为从右侧边界到该点的距离
ke = 10;                              % 衰减系数 (m^-1)，可改成 5/10/25/50 对比
y_max = max(y(:));
L = y_max - y;                        % 传播距离

% 经过衰减后的单色辐射
I = Eb .* exp(-ke .* L);

% ===================== 5. 显示辐射图（非插值） =====================
figure('Color','white');

% 用 imagesc / surf+view(2) 都行，这里用 imagesc 最直接
imagesc(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz), I);
set(gca, 'YDir', 'normal');

% 不用 jet 温度图，改成亮度图，更接近论文
colormap(gray(256));

cb = colorbar;
cb.Position = [0.73 0.15 0.03 0.7];
cb.Label.String = '单色辐射强度';
cb.Label.Rotation = 0;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.9 1.05 0];
cb.Label.FontSize = 12;

% 给颜色条加边框
pos = cb.Position;
annotation('rectangle', [pos(1)-0.03 pos(2)-0.03 pos(3)+0.1 pos(4)+0.1], ...
          'Color', 'k', 'LineWidth', 0.7);

axis tight;
axis image;
grid off;
pbaspect([1 5 1]);

xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');
title(['辐射传输成像结果, k_e = ', num2str(ke), ' m^{-1}']);

% 固定显示范围
clim([min(I(:)) max(I(:))]);
