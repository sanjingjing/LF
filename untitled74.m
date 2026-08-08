clear all;
clc;
close all;

%% ===================== 添加路径 =====================
addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\火焰频域重聚焦'));

%% ===================== 1. 基本参数 =====================
Ny = 500;
Nz = 500;

% 用于显示的物理坐标
x_axis = linspace(-0.04, 0.04, Ny);
z_axis = linspace(0, 0.4, Nz);

% 像素坐标网格
[Xp, Zp] = meshgrid(1:Ny, 1:Nz);

% 四个棋盘格在图像中的中心位置
cx1 = 150;  cz1 = 120;   % 左上
cx2 = 350;  cz2 = 120;   % 右上
cx3 = 150;  cz3 = 360;   % 左下
cx4 = 350;  cz4 = 360;   % 右下

% 四个棋盘格到主镜头的距离，单位：mm
d_list = [400, 600, 800, 1000];

%% ===================== 2. 构造四层黑白棋盘格 =====================
% 每个棋盘格总边长，单位：像素
checker_size = 100;

% 棋盘格数量：4 × 4
num_square = 4;

% 生成四个棋盘格的白格掩膜
mask1 = checkerboard_mask( ...
    Xp, Zp, cx1, cz1, checker_size, num_square);

mask2 = checkerboard_mask( ...
    Xp, Zp, cx2, cz2, checker_size, num_square);

mask3 = checkerboard_mask( ...
    Xp, Zp, cx3, cz3, checker_size, num_square);

mask4 = checkerboard_mask( ...
    Xp, Zp, cx4, cz4, checker_size, num_square);

% 直接构造灰度目标
% 白色棋格为255，黑色棋格及背景为0
obj1_gray = zeros(Nz, Ny, 'uint8');
obj2_gray = zeros(Nz, Ny, 'uint8');
obj3_gray = zeros(Nz, Ny, 'uint8');
obj4_gray = zeros(Nz, Ny, 'uint8');

obj1_gray(mask1) = 255;
obj2_gray(mask2) = 255;
obj3_gray(mask3) = 255;
obj4_gray(mask4) = 255;

%% ===================== 3. 显示四层棋盘格 =====================
figure('Color', 'white');

subplot(2,2,1);
imshow(obj1_gray, []);
title(['第1层棋盘格，d = ', num2str(d_list(1)), ' mm']);

subplot(2,2,2);
imshow(obj2_gray, []);
title(['第2层棋盘格，d = ', num2str(d_list(2)), ' mm']);

subplot(2,2,3);
imshow(obj3_gray, []);
title(['第3层棋盘格，d = ', num2str(d_list(3)), ' mm']);

subplot(2,2,4);
imshow(obj4_gray, []);
title(['第4层棋盘格，d = ', num2str(d_list(4)), ' mm']);

%% ===================== 4. 显示四层棋盘格总图 =====================
% 使用double相加，避免uint8运算中的饱和问题
obj_total_double = double(obj1_gray) + ...
                   double(obj2_gray) + ...
                   double(obj3_gray) + ...
                   double(obj4_gray);

% 限制灰度最大值为255
obj_total_gray = uint8(min(obj_total_double, 255));

figure('Color', 'white');
imagesc(x_axis, z_axis, obj_total_gray);
set(gca, 'YDir', 'normal');
colormap(gray(256));
colorbar;
axis image;
axis tight;
grid off;

xlabel('横向坐标 x (m)');
ylabel('纵向坐标 z (m)');
title('四层不同深度棋盘格总图');

%% ===================== 5. 构造RGB输入 =====================
% 当前目标是灰度图，因此三个颜色通道完全相同
obj1 = repmat(obj1_gray, [1, 1, 3]);
obj2 = repmat(obj2_gray, [1, 1, 3]);
obj3 = repmat(obj3_gray, [1, 1, 3]);
obj4 = repmat(obj4_gray, [1, 1, 3]);

%% ===================== 6. 光场相机参数 =====================
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

if mod(sen_N_total, 2) == 0
    sen_N_total = sen_N_total + 1;
end

% 目标区域的物理尺寸，单位：mm
obj_w = 0.12 * 1000;
obj_h = 0.12 * 1000;

fprintf('主镜头直径 D = %.4f mm\n', D);
fprintf('微透镜直径 lens_d = %.4f mm\n', lens_d);
fprintf('单方向微透镜数量 micr_N = %d\n', micr_N);
fprintf('传感器采样数 sen_N = %d\n', sen_N);

%% ===================== 7. 四个深度分别传播 =====================
% 四个对象都是灰度图，因此只需要计算一次，不需要重复计算RGB三通道
obje1 = double(obj1_gray);
obje2 = double(obj2_gray);
obje3 = double(obj3_gray);
obje4 = double(obj4_gray);

tic;

fprintf('\n正在计算第1层棋盘格，d = %d mm...\n', d_list(1));
im1 = LF_sim( ...
    obje1, d_list(1), obj_w, obj_h, ...
    D, N_line, F, v, lens_d, lens_f, sen_N);

fprintf('正在计算第2层棋盘格，d = %d mm...\n', d_list(2));
im2 = LF_sim( ...
    obje2, d_list(2), obj_w, obj_h, ...
    D, N_line, F, v, lens_d, lens_f, sen_N);

fprintf('正在计算第3层棋盘格，d = %d mm...\n', d_list(3));
im3 = LF_sim( ...
    obje3, d_list(3), obj_w, obj_h, ...
    D, N_line, F, v, lens_d, lens_f, sen_N);

fprintf('正在计算第4层棋盘格，d = %d mm...\n', d_list(4));
im4 = LF_sim( ...
    obje4, d_list(4), obj_w, obj_h, ...
    D, N_line, F, v, lens_d, lens_f, sen_N);

% 四个不同深度的光场直接叠加
im_sum = im1 + im2 + im3 + im4;

t = toc;

fprintf('\n四层不同深度棋盘格光场传播时间：%.3f s\n', t);

%% ===================== 8. 构造RGB光场 =====================
% 灰度目标三个通道相同
im_total_R = im_sum;
im_total_G = im_sum;
im_total_B = im_sum;

% 保存分通道数据
for k = 1:3
    save( ...
        sprintf('./dataRGB/im_dmix_v_%d_Nline_%d_%d.mat', ...
        v, N_line, k), ...
        'im_sum');
end

%% ===================== 9. 组装5D光场 =====================
if ~isequal(size(im_total_R), ...
             size(im_total_G), ...
             size(im_total_B))

    error('三个通道尺寸不一致，无法组装5D光场');
end

% 第5维表示RGB颜色通道
LF = cat(5, im_total_R, im_total_G, im_total_B);

fprintf('\n5D光场LF的尺寸为：\n');
disp(size(LF));

fprintf('LF的维数 ndims(LF) = %d\n', ndims(LF));

save( ...
    sprintf('./dataRGB/LF_5D_depthmix_v_%d_Nline_%d.mat', ...
    v, N_line), ...
    'LF', '-v7.3');

fprintf('5D光场LF保存完成！\n');

%% ===================== 10. 显示LF_sim输入总图 =====================
obj_total_rgb = repmat(obj_total_gray, [1, 1, 3]);

figure('Color', 'white');
imshow(obj_total_rgb, []);
title('输入LF\_sim的四层不同深度棋盘格总图');

%% ===================== 11. 保存兼容格式 =====================
% reshape4to2_im和sum_4D_im通过文件名读取数据
% 因此将叠加后的光场保存为原函数需要的格式

d_compat = 600;

im = im_total_R;
save( ...
    sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_compat, v, N_line, 1), ...
    'im');

im = im_total_G;
save( ...
    sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_compat, v, N_line, 2), ...
    'im');

im = im_total_B;
save( ...
    sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_compat, v, N_line, 3), ...
    'im');

%% ===================== 12. 显示微透镜阵列传感器图像 =====================
im_max = [
    max(im_total_R(:));
    max(im_total_G(:));
    max(im_total_B(:))
];

im_RGB = reshape4to2_im( ...
    d_compat, v, N_line, micr_N, sen_N, im_max);

figure('Color', 'white');
imshow(im_RGB, []);
title('微透镜阵列传感器图像（四层棋盘格叠加）');

imwrite( ...
    im_RGB, ...
    sprintf('./dataRGB/im_checkerboard_depthmix_v_%d_Nline_%d_RGB.jpg', ...
    v, N_line));

%% ===================== 13. 原始光场求和重构 =====================
im_micr = sum_4D_im( ...
    d_compat, v, N_line, micr_N);

max_im_micr = max(im_micr(:));

if max_im_micr > 0
    im_micr = uint8( ...
        double(im_micr) / double(max_im_micr) * 255);
else
    im_micr = uint8(im_micr);
end

figure('Color', 'white');
imshow(im_micr, []);
title('四层棋盘格原始数据求和重构');

fprintf('\n已画出原始数据求和重构图像！\n');
fprintf('仿真完成！\n');

disp('The end!');

%% ===================== 本地函数：生成黑白棋盘格 =====================
function mask = checkerboard_mask( ...
    Xp, Zp, cx, cz, checker_size, num_square)

    % 每个小棋格的边长
    cell_size = checker_size / num_square;

    % 棋盘格左上角坐标
    x0 = cx - checker_size / 2;
    z0 = cz - checker_size / 2;

    % 相对于棋盘格左上角的位置
    x_rel = Xp - x0;
    z_rel = Zp - z0;

    % 限制棋盘格所在区域
    inside = ...
        x_rel >= 0 & x_rel < checker_size & ...
        z_rel >= 0 & z_rel < checker_size;

    % 判断当前像素属于第几个小棋格
    ix = floor(x_rel / cell_size);
    iz = floor(z_rel / cell_size);

    % 偶数位置作为白色棋格
    white_pattern = mod(ix + iz, 2) == 0;

    % 只保留棋盘范围内的白色棋格
    mask = inside & white_pattern;
end