clear all;
clc;
close all;

hsize = 37;
num_iter = 250;

n = 512;

Q1 = imresize(im2double(rgb2gray(imread('三角形.png'))), [512 512]);
Q2 = imresize(im2double(rgb2gray(imread('正方形.png'))), [512 512]);
Q3 = imresize(im2double(rgb2gray(imread('圆形.png'))), [512 512]);

I = Q1 + Q2 + Q3;
I = mat2gray(I);

figure;
imshow(I, []);
title('叠加后的原始图像');

%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.36);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 7.55);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 12.97);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.37);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.29);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.44);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 12.72);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 7.55);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.27);

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
field_rl = field_rl ./ max(field_rl(:));

figure;
imshow(field_rl, []);
title('RL去卷积整体结果');

%% Split reconstructed layers

rec1 = field_rl(1:512, :);
rec2 = field_rl(513:1024, :);
rec3 = field_rl(1025:1536, :);

figure;
subplot(1,3,1); imshow(rec1, []); title('恢复层 x1');
subplot(1,3,2); imshow(rec2, []); title('恢复层 x2');
subplot(1,3,3); imshow(rec3, []); title('恢复层 x3');