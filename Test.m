clear;
clc;
close all;

%% ===================== 1. 添加路径 =====================
addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\火焰频域重聚焦'));
%% ===================== 2. 读取 RGB 原图 =====================
I1 = imread('LenaRGB.tif');
I2 = imread('BaboonRGB.tif');
I3 = imread('PeppersRGB.tif');
figure;
subplot(1,3,1),imshow(I1);
subplot(1,3,2),imshow(I2);
subplot(1,3,3),imshow(I3);
d_list = [400,600, 800];% 三层到主镜头距离（mm）
% 转成 double，但保留原来的 0～255 强度范围
obj1 = double(I1);
obj2 = double(I2);
obj3 = double(I3);
%% ===================== 3. 光场相机参数 =====================
v = 16;
N_line = 81;
D = 4;
F = 16;
lens_d = 0.01;
sen_N = 20;
lens_f = lens_d * F / D;
lens_v = lens_f + v;
micr_N = ceil(D / lens_d);
sen_N_total = micr_N * sen_N;
sen_d = lens_d / sen_N;

if mod(sen_N_total,2) == 0
    sen_N_total = sen_N_total + 1;
end

obj_w = 0.12 * 1000;   % mm
obj_h = 0.12 * 1000;   % mm

%% ===================== 9. 分层传播并叠加 =====================
im_total_R = [];
im_total_G = [];
im_total_B = [];

tic
for k = 1:3
    fprintf('正在计算第 %d 个通道...\n', k);

    obje1 = double(obj1(:,:,k));
    obje2 = double(obj2(:,:,k));
    obje3 = double(obj3(:,:,k));
    obje4 = double(obj4(:,:,k));

    im1 = LF_sim(obje1, d_list(1), obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im2 = LF_sim(obje2, d_list(2), obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im3 = LF_sim(obje3, d_list(3), obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im4 = LF_sim(obje4, d_list(4), obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);

    im_sum = im1 + im2 + im3 + im4;

    save(sprintf('./dataRGB/im_dmix_v_%d_Nline_%d_%d.mat', v, N_line, k), 'im_sum');

    if k == 1
        im_total_R = im_sum;
    elseif k == 2
        im_total_G = im_sum;
    else
        im_total_B = im_sum;
    end
end
t = toc;
fprintf('\n模拟四层不同深度光场传输时间 t = %.3f s\n', t);

%% ===================== 10. 组装 5D 光场 =====================
if ~isequal(size(im_total_R), size(im_total_G), size(im_total_B))
    error('三个通道尺寸不一致，无法组装 5D 光场');
end

LF = cat(5, im_total_R, im_total_G, im_total_B);

fprintf('\n5D 光场 LF 的维度为: ');
disp(size(LF));
fprintf('LF 的维数 ndims(LF) = %d\n', ndims(LF));

save(sprintf('./dataRGB/LF_5D_depthmix_v_%d_Nline_%d.mat', v, N_line), 'LF', '-v7.3');
fprintf('5D 光场 LF 已保存完成！\n');

%% ===================== 11. 显示输入场景（总图） =====================
obj_total = double(obj1_gray) + double(obj2_gray) + double(obj3_gray) + double(obj4_gray);
max_val = double(max(obj_total(:)));

if max_val > 0
    obj_total = uint8(obj_total / max_val * 255);
else
    obj_total = uint8(obj_total);
end
obj_total_rgb = repmat(obj_total, [1,1,3]);

figure;
imshow(obj_total_rgb, []);
title('输入 LF\_sim 的四层不同深度图形总图');

%% ===================== 12. 显示微透镜阵列传感器图像 =====================
% 这里需要 im_max，因此从叠加后的三个通道计算
im_max = [max(im_total_R(:)); max(im_total_G(:)); max(im_total_B(:))];

% 为了兼容 reshape4to2_im，临时保存成原函数读取格式
im = im_total_R; save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', 600, v, N_line, 1), 'im');
im = im_total_G; save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', 600, v, N_line, 2), 'im');
im = im_total_B; save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', 600, v, N_line, 3), 'im');

im_RGB = reshape4to2_im(600, v, N_line, micr_N, sen_N, im_max);

figure;
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像（四层不同深度叠加）');

imwrite(im_RGB, sprintf('./dataRGB/im_depthmix_v_%d_Nline_%d_RGB.jpg', v, N_line));

fprintf('\n仿真完成：已考虑四个图形距离主镜头不同。\n');
d=600;
%% 原始数据求和重构后的图像
im_micr=sum_4D_im(d,v,N_line,micr_N);
max_im_micr=max(max(max(im_micr)));
im_micr=uint8(im_micr/max_im_micr*255);
figure
imshow(im_micr,[]);title('原始数据求和重构：');
fprintf('\n已画出原始数据重构图像!\n');
disp('The end!');
