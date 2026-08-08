clear all;
clc;
close all;

hsize = 37;
num_iter = 50;

n = 512;

Q1 = imresize(im2double(rgb2gray(imread('flame1.bmp'))), [512 512]);
Q2 = imresize(im2double(rgb2gray(imread('flame100.bmp'))), [512 512]);
Q3 = imresize(im2double(rgb2gray(imread('flame191.bmp'))), [512 512]);



%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 2.36);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 5.55);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 7.97);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 5.37);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 2.29);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 5.44);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 7.72);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 5.55);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator(size(Q1), hsize, 2.27);



%% Stack measured images b = [b1; b2; b3]

image = zeros(1536, 512);
image(1:512, :)       = Q1;
image(513:1024, :)    = Q2;
image(1025:1536, :)   = Q3;

image = max(real(image), 0);

figure;
imshow(mat2gray(image));
title('观测图像');

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
subplot(1,3,1); imshow(Q1, []); title('重聚焦 x1');
subplot(1,3,2); imshow(Q2, []); title('重聚焦 x2');
subplot(1,3,3); imshow(Q3, []); title('重聚焦 x3');
figure;
subplot(1,3,1); imshow(rec1, []); title('恢复层 x1');
subplot(1,3,2); imshow(rec2, []); title('恢复层 x2');
subplot(1,3,3); imshow(rec3, []); title('恢复层 x3');