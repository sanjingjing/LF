clear ;
clc;
close all;

hsize = 37;
num_iter = 50;
image_size=512;
imageFolder = fullfile(pwd, '重聚焦图像');

%% ===================== 1. 火焰辐射模型 =====================
Q1 = imread(fullfile( imageFolder, '043_alpha_0.42.bmp'));

Q2 = imread(fullfile(imageFolder, '044_alpha_0.43.bmp'));
Q3 = imread(fullfile(imageFolder,  '045_alpha_0.44.bmp'));

%% 转换为double类型
Q1 = im2double(im2gray(Q1));
Q2 = im2double(im2gray(Q2));
Q3 =im2double(im2gray(Q3));
Q1 = imresize(Q1,[image_size,image_size]);
Q2 =imresize(Q2,[image_size,image_size]);
Q3 =imresize(Q3,[image_size,image_size]);


%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 3.37);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 6.81);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 9.73);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 6.36);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 2.17);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 6.59);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 9.05);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 6.19);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator([image_size,image_size], hsize, 2.61);



%% Stack measured images b = [b1; b2; b3]

image = zeros(1536, 512);
image(1:512, :)       = Q1;
image(513:1024, :)    = Q2;
image(1025:1536, :)   = Q3;

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
rec1 = field_rl(1:512, :);
rec2 = field_rl(513:1024, :);
rec3 = field_rl(1025:1536, :);
figure;
subplot(2,3,1);imshow(Q1, []);title('重聚焦层1');
subplot(2,3,2);imshow(Q2, []);title('重聚焦层2');
subplot(2,3,3);imshow(Q3, []);title('重聚焦层3');
subplot(2,3,4);imshow(rec1, []);title('恢复层 x1');
subplot(2,3,5);imshow(rec2, []);title('恢复层 x2');
subplot(2,3,6);imshow(rec3, []);title('恢复层 x3');

