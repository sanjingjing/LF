clear all;
clc; 
close all;
% ===================== 1. 定义模型常数 =====================
t1 = 1200;   % 温度增量 (K)
t2 = 900;    % 环境温度 (K)
m = 3;       % 衰减系数
n = 0.9;       % 形状参数
R=0.04;
Z=0.4;
x=0;
[y,z] = meshgrid(linspace(-0.4, 0.4, 500),linspace(0, 0.4, 500));
term = m*( x.^2/R.^2 +y.^2/R.^2+ z.^2/Z.^2 )-n; % 指数内的核心项

T = t1 * exp(-1* term.^2) + t2;        % 温度场计算
figure('Color','white');
surf(y,z,T);axis tight;
axis normal;  % MATLAB 根据数据自动缩放，不强制像素比例
colormap(jet(11));shading interp;         % 平滑颜色

