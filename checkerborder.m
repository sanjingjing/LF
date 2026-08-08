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
Ny = 500;                     % 图像宽度，像素
Nz = 500;                     % 图像高度，像素

% 灰度像素值
pixel_black = uint8(0);       % 黑色和背景
pixel_white = uint8(255);     % 白色棋盘格

% 像素坐标网格
[Xp, Zp] = meshgrid(1:Ny, 1:Nz);

% 四个棋盘格在图像中的中心位置
cx1 = 150;
cz1 = 120;                    % 第1层：左上

cx2 = 350;
cz2 = 120;                    % 第2层：右上

cx3 = 150;
cz3 = 360;                    % 第3层：左下

cx4 = 350;
cz4 = 360;                    % 第4层：右下

% 四层目标到主镜头的距离，单位：mm
d_list = [400, 600, 800, 1000];

%% ===================== 2. 构造四层棋盘格像素图 =====================
% 每层直接使用0～255灰度像素，不再使用温度、辐射和衰减模型
P1 = zeros(Nz, Ny, 'uint8');
P2 = zeros(Nz, Ny, 'uint8');
P3 = zeros(Nz, Ny, 'uint8');
P4 = zeros(Nz, Ny, 'uint8');

% 棋盘格参数
checker_size = 100;           % 棋盘格总边长，单位：像素
num_square = 4;               % 4×4棋盘格

% 棋盘格掩膜生成函数
make_checker = @(cx, cz) ...
    checkerboard_mask(Xp, Zp, cx, cz, checker_size, num_square);

% 生成四层棋盘格
mask1 = make_checker(cx1, cz1);
mask2 = make_checker(cx2, cz2);
mask3 = make_checker(cx3, cz3);
mask4 = make_checker(cx4, cz4);

% 黑格和背景为0，白格为255
P1(:) = pixel_black;
P2(:) = pixel_black;
P3(:) = pixel_black;
P4(:) = pixel_black;

P1(mask1) = pixel_white;
P2(mask2) = pixel_white;
P3(mask3) = pixel_white;
P4(mask4) = pixel_white;

% LF_sim输入图像
obj1_gray = P1;
obj2_gray = P2;
obj3_gray = P3;
obj4_gray = P4;

% 转为三通道灰度图，兼容原来的RGB计算过程
obj1 = repmat(obj1_gray, [1, 1, 3]);
obj2 = repmat(obj2_gray, [1, 1, 3]);
obj3 = repmat(obj3_gray, [1, 1, 3]);
obj4 = repmat(obj4_gray, [1, 1, 3]);

%% ===================== 3. 显示各层棋盘格 =====================
figure('Color', 'white');

subplot(2, 2, 1);
imshow(obj1_gray, [0, 255]);
title(['第1层棋盘格，d = ', num2str(d_list(1)), ' mm']);

subplot(2, 2, 2);
imshow(obj2_gray, [0, 255]);
title(['第2层棋盘格，d = ', num2str(d_list(2)), ' mm']);

subplot(2, 2, 3);
imshow(obj3_gray, [0, 255]);
title(['第3层棋盘格，d = ', num2str(d_list(3)), ' mm']);

subplot(2, 2, 4);
imshow(obj4_gray, [0, 255]);
title(['第4层棋盘格，d = ', num2str(d_list(4)), ' mm']);

%% ===================== 4. 显示四层像素总图 =====================
% 使用最大值进行合成，避免像素叠加后超过255
pixel_total = max( ...
    cat(3, obj1_gray, obj2_gray, obj3_gray, obj4_gray), ...
    [], 3);

figure('Color', 'white');
imshow(pixel_total, [0, 255]);
axis image;
xlabel('x / pixel');
ylabel('z / pixel');
title('四层不同深度棋盘格像素总图');

%% ===================== 5. 光场相机参数 =====================
v = 16;                       % 主透镜到微透镜阵列的距离（主透镜像距）/mm
N_line = 81;                  % 光线采样数

D = 4;                        % 主镜头孔径/mm
F = 16;                       % 主镜头焦距/mm
lens_d = 0.01;                % 微透镜直径/mm
sen_N = 20;                   % 单个微透镜对应的像素数

lens_f = lens_d * F / D;      % 微透镜焦距，即微透镜阵列到传感器的距离/mm
lens_v = lens_f + v;           % 主透镜到图像传感器的总距离/mm

micr_N = ceil(D / lens_d);    % 单方向微透镜数量
sen_N_total = micr_N * sen_N; % 单方向传感器总像素数
sen_d = lens_d / sen_N;       % 图像传感器像素尺寸/mm

if mod(sen_N_total, 2) == 0
    sen_N_total = sen_N_total + 1;
end

% 目标实际尺寸，单位：mm
obj_w = 0.12 * 1000;
obj_h = 0.12 * 1000;

%% ===================== 6. 分层传播并叠加 =====================
im_total_R = [];
im_total_G = [];
im_total_B = [];

tic;

for k = 1:3

    fprintf('正在计算第 %d 个通道...\n', k);

    % 当前颜色通道的四层像素图
    obje1 = double(obj1(:, :, k));
    obje2 = double(obj2(:, :, k));
    obje3 = double(obj3(:, :, k));
    obje4 = double(obj4(:, :, k));

    % 四个不同深度分别进行光场传播
    im1 = LF_sim( ...
        obje1, d_list(1), obj_w, obj_h, ...
        D, N_line, F, v, lens_d, lens_f, sen_N);

    im2 = LF_sim( ...
        obje2, d_list(2), obj_w, obj_h, ...
        D, N_line, F, v, lens_d, lens_f, sen_N);

    im3 = LF_sim( ...
        obje3, d_list(3), obj_w, obj_h, ...
        D, N_line, F, v, lens_d, lens_f, sen_N);

    im4 = LF_sim( ...
        obje4, d_list(4), obj_w, obj_h, ...
        D, N_line, F, v, lens_d, lens_f, sen_N);

    % 四层光场叠加
    im_sum = im1 + im2 + im3 + im4;

    save( ...
        sprintf( ...
        './dataRGB/im_dmix_v_%d_Nline_%d_%d.mat', ...
        v, N_line, k), ...
        'im_sum');

    % 保存三个通道
    if k == 1
        im_total_R = im_sum;

    elseif k == 2
        im_total_G = im_sum;

    else
        im_total_B = im_sum;
    end
end

t = toc;

fprintf( ...
    '\n模拟四层不同深度光场传输时间 t = %.3f s\n', ...
    t);

%% ===================== 7. 组装5D光场 =====================
if ~isequal( ...
        size(im_total_R), ...
        size(im_total_G), ...
        size(im_total_B))

    error('三个通道尺寸不一致，无法组装5D光场');
end

LF = cat(5, im_total_R, im_total_G, im_total_B);

fprintf('\n5D光场LF的维度为：');
disp(size(LF));

fprintf('LF的维数 ndims(LF) = %d\n', ndims(LF));

save( ...
    sprintf( ...
    './dataRGB/LF_5D_depthmix_v_%d_Nline_%d.mat', ...
    v, N_line), ...
    'LF', ...
    '-v7.3');

fprintf('5D光场LF已保存完成！\n');

%% ===================== 8. 显示LF_sim输入场景 =====================
obj_total_rgb = repmat(pixel_total, [1, 1, 3]);

figure('Color', 'white');
imshow(obj_total_rgb, [0, 255]);
title('输入LF\_sim的四层不同深度棋盘格');

%% ===================== 9. 显示微透镜阵列传感器图像 =====================
% 计算三个通道的最大值
im_max = [
    max(im_total_R(:));
    max(im_total_G(:));
    max(im_total_B(:))
];

% reshape4to2_im和sum_4D_im使用的参考距离标记
d_ref = 600;

% 按原函数要求保存为对应的文件名称
im = im_total_R;

save( ...
    sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_ref, v, N_line, 1), ...
    'im');

im = im_total_G;

save( ...
    sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_ref, v, N_line, 2), ...
    'im');

im = im_total_B;

save( ...
    sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_ref, v, N_line, 3), ...
    'im');

% 将4D微透镜数据转换成2D传感器图像
im_RGB = reshape4to2_im( ...
    d_ref, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure('Color', 'white');
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像');

imwrite( ...
    im_RGB, ...
    sprintf( ...
    './dataRGB/im_depthmix_v_%d_Nline_%d_RGB.jpg', ...
    v, N_line));

fprintf( ...
    '\n仿真完成：输入采用0～255灰度像素，' ...
    );
fprintf('已考虑四层目标距离主镜头不同。\n');

%% ===================== 10. 原始数据求和重构 =====================
im_micr = sum_4D_im( ...
    d_ref, ...
    v, ...
    N_line, ...
    micr_N);

max_im_micr = max(im_micr(:));

if max_im_micr > 0
    im_micr = uint8( ...
        double(im_micr) / double(max_im_micr) * 255);
else
    im_micr = uint8(im_micr);
end

figure('Color', 'white');
imshow(im_micr, []);
title('原始数据求和重构');

fprintf('\n已画出原始数据重构图像！\n');
disp('The end!');

%% ===================== 本地函数：生成棋盘格 =====================
function mask = checkerboard_mask( ...
    Xp, Zp, cx, cz, checker_size, num_square)

    % 单个方格的边长，单位：像素
    cell_size = checker_size / num_square;

    % 棋盘格左上角位置
    x0 = cx - checker_size / 2;
    z0 = cz - checker_size / 2;

    % 相对于棋盘格左上角的像素坐标
    x_rel = Xp - x0;
    z_rel = Zp - z0;

    % 判断像素是否位于棋盘格范围内
    inside = ...
        x_rel >= 0 & ...
        x_rel < checker_size & ...
        z_rel >= 0 & ...
        z_rel < checker_size;

    % 计算当前像素属于第几个方格
    ix = floor(x_rel / cell_size);
    iz = floor(z_rel / cell_size);

    % 黑白交替，白格为true
    pattern = mod(ix + iz, 2) == 0;

    % 只保留棋盘格范围内的白格
    mask = inside & pattern;
end