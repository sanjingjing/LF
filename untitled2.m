clear all;
clc;
close all;
hsize=37;
max_photons = 100;
num_iter = 30;
pixel_size = 50;
spacing_px = 4;
n = 512;
lambda = 510;
numerical_aperture = 1;
background_level = 0;
left_bg = 0;
mid_bg = 0.05;
right_bg = 0.25;

SAVE_DISPLAY = 1;
USE_GPU = 0;

Q1 = imresize(im2double(rgb2gray(imread('三角形.png'))), [512 512]);
Q2 = imresize(im2double(rgb2gray(imread('正方形.png'))), [512 512]);
Q3 = imresize(im2double(rgb2gray(imread('圆形.png'))), [512 512]);

I = Q1 + Q2 + Q3;
I = mat2gray(I);  % 归一化到 [0,1]，防止亮度过低显示空白
figure;
imshow(I); 
title('叠加后的原始图像');
%% Create PSF
otf1 = paraxial_otf(n, lambda, numerical_aperture, pixel_size);


 [Op11, ~, fftFilter11] = createGaussianBlurringOperator(size(I),hsize,2.36);
[Op12, ~, fftFilter12] = createGaussianBlurringOperator(size(I),hsize,3.55);
[Op13, ~, fftFilter13] = createGaussianBlurringOperator(size(I),hsize,5.97);
[Op21, ~, fftFilter21] = createGaussianBlurringOperator(size(I),hsize,3.37);
[Op22, ~, fftFilter22] = createGaussianBlurringOperator(size(I),hsize,2.29);
[Op23, ~, fftFilter23] = createGaussianBlurringOperator(size(I),hsize,3.44);
[Op31, ~, fftFilter31] = createGaussianBlurringOperator(size(I),hsize,5.72);
[Op32, ~, fftFilter32] = createGaussianBlurringOperator(size(I),hsize,3.25);
[Op33, ~, fftFilter33] = createGaussianBlurringOperator(size(I),hsize,2.27);


I1=Op11(Q1)+Op12(Q2)+Op13(Q3);
I2=Op21(Q1)+Op22(Q2)+Op23(Q3);
I3=Op31(Q1)+Op32(Q2)+Op33(Q3);
figure;
I1 = mat2gray(I1);  
I2 = mat2gray(I2); 
I3 = mat2gray(I3); 
subplot(1,3,1); imshow(I1, []); title('叠加图像1');
subplot(1,3,2); imshow(I2, []); title('叠加图像2');
subplot(1,3,3); imshow(I3, []); title('叠加图像3');


image = zeros(1536,512);
image(1:512, :)       = I1;
image(513:1024, :)    = I2;
image(1025:1536, :)   = I3;
figure;
imshow(image);

fftFilter = zeros(1536,1536);

fftFilter(1:512,       1:512)       = fftFilter11;
fftFilter(1:512,       513:1024)    = fftFilter12;
fftFilter(1:512,       1025:1536)   = fftFilter13;

fftFilter(513:1024,    1:512)       = fftFilter21;
fftFilter(513:1024,    513:1024)    = fftFilter22;
fftFilter(513:1024,    1025:1536)   = fftFilter23;

fftFilter(1025:1536,   1:512)       = fftFilter31;
fftFilter(1025:1536,   513:1024)    = fftFilter32;
fftFilter(1025:1536,   1025:1536)   = fftFilter33;

field_rl = richardson_lucy(image, fftFilter, num_iter, 1);
 field_rl = field_rl ./ max(field_rl(:));

 
 figure;

imshow(field_rl, []); title('RL去卷积图像');









