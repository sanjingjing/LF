clear all;
clc;
close all;

hsize = 37;
num_iter = 10;
image_size= 512;

%% ===================== 1. 火焰辐射模型 =====================
t1 = 1200;     % 温度增量 (K)
t2 = 900;      % 环境温度 (K)
m  = 3;        % 形状参数
n  = 0.9;
R  = 0.040;     % m
Z  = 0.4;      % m

x1 = 0;  % 固定切面位置





Ny = 500;
Nz = 500;
[y,z] = meshgrid(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz));
% 温度场
term1 = m * (x1.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;

T1 = t1 * exp(-(term1.^2)) + t2;


% 单色辐射模型
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K
Eb1 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T1)) - 1);


% Beer 衰减
ke = 8;             % m^-1
x_max = 0.025;
L1 = x_max - x1;


Q1 = Eb1 .* exp(-ke .* L1);



Q1 = flipud(imresize(im2double(Q1), [512 512]));



I = Q1 ;
I = mat2gray(I);

figure;
imshow(I, []);
title('原始图像');


