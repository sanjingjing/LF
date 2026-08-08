clear all;
clc;
close all;


addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\频域重聚焦精简'));
%% ===================== 1. 火焰辐射模型 =====================
t1 = 1200;     % 温度增量 (K)
t2 = 900;      % 环境温度 (K)
m  = 3;        % 形状参数
n  = 0.9;
R  = 0.04;     % m
Z  = 0.4;      % m
x  = 0.02;     % 固定切面位置

Ny = 500;
Nz = 500;
[y,z] = meshgrid(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz));

% 温度场
term = m * (x.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
T = t1 * exp(-(term.^2)) + t2;

% 单色辐射模型
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K
Eb = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T)) - 1);

% Beer 衰减
ke = 30;             % m^-1
y_max = max(y(:));
L = y_max - y;
I = Eb .* exp(-ke .* L);

%% ===================== 2. 显示火焰辐射图 =====================
figure('Color','white');
imagesc(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz), I);
set(gca, 'YDir', 'normal');
colormap(gray(256));

cb = colorbar;
cb.Position = [0.73 0.15 0.03 0.7];
cb.Label.String = '单色辐射强度';
cb.Label.Rotation = 0;
cb.Label.Units = 'normalized';
cb.Label.Position = [0.9 1.05 0];
cb.Label.FontSize = 12;

axis tight;
axis image;
grid off;
pbaspect([1 5 1]);
xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');
title(['辐射传输成像结果, k_e = ', num2str(ke), ' m^{-1}']);
clim([min(I(:)) max(I(:))]);

%% ===================== 3. 把火焰辐射图转成 LF_sim 输入 =====================
I_norm = I - min(I(:));
if max(I_norm(:)) > 0
    I_norm = I_norm / max(I_norm(:));
end

I_norm = I_norm .^ 0.7;   % gamma 调整

obj_gray = uint8(I_norm * 255);
obj = repmat(obj_gray, [1, 1, 3]);   % RGB 三通道

object_num = 1;
d = 600;   % mm，火焰到主透镜距离

%% ===================== 4. 光场相机参数 =====================
v = 16;          % 像距
N_line = 81;     % 主透镜离散化
d_m = 18;

D = 4;           % 主透镜直径
F = 16;          % 焦距
lens_d = 0.02;   % 微透镜直径
sen_N = 20;      % 每个微透镜后传感器像素数

lens_f = lens_d * F / D;
lens_v = lens_f + v;
micr_N = ceil(D / lens_d);
sen_N_total = micr_N * sen_N;
sen_d = lens_d / sen_N;

if mod(sen_N_total,2) == 0
    sen_N_total = sen_N_total + 1;
end

im_max = zeros(3,1);

%% ===================== 5. 火焰辐射图送入 LF_sim =====================
tic
for k = 1:3
    fprintf('正在计算第 %d 个通道...\n', k);
    obje = double(obj(:,:,k));
   obj_w = 0.12 * 1000;   % 横向物理宽度，0.08 m = 80 mm
obj_h = 0.60 * 1000;   % 纵向物理高度，0.40 m = 400 mm

    im = LF_sim(obje, d, obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im_max(k) = max(im(:));

    save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, k), 'im');
end
t = toc;
fprintf('\n模拟光场传输时间为：t = %.3f s\n', t);
%% ===================== 5.1 组装 5D 光场数据 =====================
LF_R = load(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, 1));
LF_G = load(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, 2));
LF_B = load(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, 3));

LF_R = LF_R.im;
LF_G = LF_G.im;
LF_B = LF_B.im;

fprintf('\nR通道维度: ');
disp(size(LF_R));
fprintf('G通道维度: ');
disp(size(LF_G));
fprintf('B通道维度: ');
disp(size(LF_B));

% 先检查三个通道尺寸是否一致
if ~isequal(size(LF_R), size(LF_G), size(LF_B))
    error('三个通道的光场尺寸不一致，不能合并成5D光场。');
end

% 组装为 5D 光场
LF = cat(5, LF_R, LF_G, LF_B);

fprintf('\n5D光场 LF 的维度为: ');
disp(size(LF));
fprintf('LF 的维数 ndims(LF) = %d\n', ndims(LF));

save(sprintf('./dataRGB/LF_5D_d_%d_v_%d_Nline_%d.mat', d, v, N_line), 'LF', '-v7.3');

fprintf('5D光场数据已保存完成！\n');

%% ===================== 6. 显示输入场景 =====================
figure;
imshow(obj, []);
title('输入 LF_sim 的火焰辐射图 obj');
fprintf('\n已画出火焰场景图像！\n');

%% ===================== 7. 显示传感器图像 =====================
im_RGB = reshape4to2_im(d, v, N_line, micr_N, sen_N, im_max);

figure;
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像');
imwrite(im_RGB, sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_RGB.jpg', d, v, N_line));
% 设置显示比例，确保纵横比与火焰图像一致
pbaspect([1 5 1]);  % 调整纵横比为 1:5，匹配火焰图像
fprintf('\n已画出传感器图像！\n');

%% ===================== 8. 原始数据求和重构 =====================
im_micr = sum_4D_im(d, v, N_line, micr_N);
max_im_micr = max(im_micr(:));
im_micr = uint8(im_micr / max_im_micr * 255);

figure;
imshow(im_micr, []);
title('原始数据求和重构');
fprintf('\n已画出原始数据重构图像!\n');

disp('The end!');
