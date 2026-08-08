clear all;
clc;
close all;
max_photons = 100;
num_iter = 100;
pixel_size = 20;
spacing_px = 4;
n = 256;
lambda = 510;
numerical_aperture = 2.4;
background_level = 0;
left_bg = 0;
mid_bg = 0.05;
right_bg = 0.25;

SAVE_DISPLAY = 1;
USE_GPU = 0;
I=im2double(imread("cameraman.tif"));

%% Create OTF
otf = paraxial_otf(n, lambda, numerical_aperture, pixel_size);







%% Simulate captured data
field_imaged = real(ifft2(fft2(I) .* otf));
field_imaged = field_imaged ./ max(field_imaged(:));

I_use= poissrnd(field_imaged * max_photons + background_level);

 
 field_rl = richardson_lucy(I_use, otf, num_iter, 1);
 field_rl = field_rl ./ max(field_rl(:));

 
 figure;
subplot(1,3,1); imshow(I, []); title('原始图像');
subplot(1,3,2); imshow(I_use, []); title('模糊图像');
subplot(1,3,3); imshow(field_rl, []); title('RL去卷积图像');








