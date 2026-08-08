clear all;
clc;
close all;

%% ===================== 1. 火焰辐射模型 =====================
t1 = 1200;     % 温度增量 (K)
t2 = 900;      % 环境温度 (K)
m  = 10;        % 形状参数
n  = 0.4;
R  = 0.025;     % m
Z  = 0.075;      % m

x1 = 0;  % 固定切面位置
x2 = 0.004; 
x3 = 0.007; 
x4 = -0.007; 



Ny = 500;
Nz = 500;
[y,z] = meshgrid(linspace(-0.02, 0.02, Ny), linspace(0, 0.075, Nz));
% 温度场
term1 = m * (x1.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term2 = m * (x2.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term3 = m * (x3.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term4 = m * (x4.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
T1 = t1 * exp(-(term1.^2)) + t2;
T2 = t1 * exp(-(term2.^2)) + t2;
T3 = t1 * exp(-(term3.^2)) + t2;
T4 = t1 * exp(-(term4.^2)) + t2;
% 单色辐射模型
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K
Eb1 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T1)) - 1);
Eb2 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T2)) - 1);
Eb3 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T3)) - 1);
Eb4 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T4)) - 1);
% Beer 衰减
ke = 30;             % m^-1
x_max = 0.025;
L1 = x_max - x1;
L2 = x_max - x2;
L3= x_max - x3;
L4 = x_max - x4;
I1 = Eb1 .* exp(-ke .* L1);
I2 = Eb2 .* exp(-ke .* L2);
I3 = Eb3 .* exp(-ke .* L3);
I4 = Eb4 .* exp(-ke .* L4);
I=I1+I2;


%% ===================== 2. 显示火焰辐射图 =====================
figure('Color','white');
imagesc(linspace(-0.02, 0.02, Ny), linspace(0, 0.075, Nz), I);
set(gca, 'YDir', 'normal');
colormap(hot(256));
cb = colorbar;
cb.Position = [0.93 0.15 0.03 0.7];
cb.Label.String = '单色辐射强度';
cb.Label.Rotation = 0;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.75 1.15 3.5];
cb.Label.FontSize = 12;
axis tight;
axis image;
grid off;
xlabel('径向坐标 y (m)');
ylabel('轴向坐标 z (m)');
title(['辐射传输成像结果, k_e = ', num2str(ke), ' m^{-1}']);
clim([min(I(:)) max(I(:))]);





