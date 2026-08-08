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
R  = 0.025;     % m
Z  = 0.075;      % m

x1 = 0;  % 固定切面位置
x2 = 0.005; 
x3 = 0.010; 




Ny = 500;
Nz = 500;
[y,z] = meshgrid(linspace(-0.02, 0.02, Ny), linspace(0, 0.075, Nz));
% 温度场
term1 = m * (x1.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term2 = m * (x2.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term3 = m * (x3.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;

T1 = t1 * exp(-(term1.^2)) + t2;
T2 = t1 * exp(-(term2.^2)) + t2;
T3 = t1 * exp(-(term3.^2)) + t2;

% 单色辐射模型
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K
Eb1 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T1)) - 1);
Eb2 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T2)) - 1);
Eb3 = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T3)) - 1);

% Beer 衰减
ke = 8;             % m^-1
x_max = 0.025;
L1 = x_max - x1;
L2 = x_max - x2;
L3= x_max - x3;

Q1 = Eb1 .* exp(-ke .* L1);
Q2 = Eb2 .* exp(-ke .* L2);
Q3 = Eb3 .* exp(-ke .* L3);


Q1 = flipud(imresize(im2double(Q1), [512 512]));
Q2 = flipud(imresize(im2double(Q2), [512 512]));
Q3 = flipud(imresize(im2double(Q3), [512 512]));

I = Q1 + Q2 + Q3;
I = mat2gray(I);

figure;
imshow(I, []);
title('叠加后的原始图像');

figure;
subplot(1,3,1); imshow(mat2gray(Q1)); title('原始图像1');
subplot(1,3,2); imshow(mat2gray(Q2)); title('原始图像2');
subplot(1,3,3); imshow(mat2gray(Q3)); title('原始图像3');

%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.29);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.81);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 10.73);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.36);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.17);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.59);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 9.05);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.19);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.61);

%% Forward imaging model

I1 = Op11(Q1) + Op12(Q2) + Op13(Q3);
I2 = Op21(Q1) + Op22(Q2) + Op23(Q3);
I3 = Op31(Q1) + Op32(Q2) + Op33(Q3);

figure;
subplot(1,3,1); imshow(mat2gray(I1)); title('叠加图像1');
subplot(1,3,2); imshow(mat2gray(I2)); title('叠加图像2');
subplot(1,3,3); imshow(mat2gray(I3)); title('叠加图像3');

%% Stack measured images b = [b1; b2; b3]

image = zeros(1536, 512);
image(1:512, :)       = I1;
image(513:1024, :)    = I2;
image(1025:1536, :)   = I3;

image = max(real(image), 0);

figure;
imshow(mat2gray(image));
title('观测叠加图像 b');

%% Define block forward operator A

A = @(x) [
    Op11(x(1:512, :)) + Op12(x(513:1024, :)) + Op13(x(1025:1536, :));
    Op21(x(1:512, :)) + Op22(x(513:1024, :)) + Op23(x(1025:1536, :));
    Op31(x(1:512, :)) + Op32(x(513:1024, :)) + Op33(x(1025:1536, :))
];

%% Define adjoint operator A_adj

A_adj = @(y) [
    Op11_adj(y(1:512, :)) + Op21_adj(y(513:1024, :)) + Op31_adj(y(1025:1536, :));
    Op12_adj(y(1:512, :)) + Op22_adj(y(513:1024, :)) + Op32_adj(y(1025:1536, :));
    Op13_adj(y(1:512, :)) + Op23_adj(y(513:1024, :)) + Op33_adj(y(1025:1536, :))
];

%% Initial estimate

init = ones(size(image)) .* mean(image(:));

%% Richardson-Lucy deconvolution

field_rl = deconvRL_operator(A, A_adj, image, num_iter, init);

field_rl = real(field_rl);
field_rl(field_rl < 0) = 0;


figure;
imshow(field_rl, []);
title('RL去卷积整体结果');

%% Split reconstructed layers
vmin = min(field_rl(:));
vmax = max(field_rl(:));
rec1 = field_rl(1:512, :);
rec2 = field_rl(513:1024, :);
rec3 = field_rl(1025:1536, :);

figure;
subplot(1,3,1); imshow(rec1, [vmin vmax]); title('恢复层 x1');
subplot(1,3,2); imshow(rec2, [vmin vmax]); title('恢复层 x2');
subplot(1,3,3); imshow(rec3, [vmin vmax]); title('恢复层 x3');