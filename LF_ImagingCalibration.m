clear;
clc;
close all;

%% ============================================================
%  添加程序路径
%% ============================================================

addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\频域重聚焦精简'));

% 检查输出文件夹
if ~exist('./dataRGB', 'dir')
    mkdir('./dataRGB');
end

%% ============================================================
%  1. 生成黑白阶跃标定板
%
%  左半部分：黑色，灰度值0
%  右半部分：白色，灰度值255
%% ============================================================

board_rows = 500;       % 标定板图像行数
board_cols = 500;       % 标定板图像列数

% 阶跃边缘所在列，默认位于标定板中心
edge_col = floor(board_cols / 2);

% 初始化黑色标定板
obj_gray = zeros(board_rows, board_cols, 'uint8');

% 右半部分设置为白色
obj_gray(:, edge_col + 1:end) = 255;

% 转换为RGB三通道
obj = repmat(obj_gray, [1, 1, 3]);

fprintf('标定板图像尺寸：\n');
disp(size(obj));

fprintf('阶跃边缘所在列：%d\n', edge_col);
fprintf('左侧灰度值：%d\n', obj_gray(1, 1));
fprintf('右侧灰度值：%d\n', obj_gray(1, end));

%% ============================================================
%  2. 显示理想黑白标定板
%% ============================================================

figure('Color', 'white');

imshow(obj_gray, []);

title('光学系统输入：黑白阶跃标定板');
xlabel('横向像素');
ylabel('纵向像素');

%% ============================================================
%  3. 显示标定板中心行的理想阶跃响应
%% ============================================================

center_row = round(board_rows / 2);
ideal_edge_profile = double(obj_gray(center_row, :)) / 255;

figure('Color', 'white');

plot(1:board_cols, ideal_edge_profile, ...
    'LineWidth', 1.5);

xline(edge_col + 0.5, '--', '阶跃边缘');

xlabel('横向像素');
ylabel('归一化灰度');
title('黑白标定板理想阶跃函数');
ylim([-0.1, 1.1]);
grid on;
box on;

%% ============================================================
%  4. 标定板物理参数
%% ============================================================

% 标定板中心到主透镜的距离
d = 603;              % mm

% 标定板实际物理尺寸
%
% 为避免原火焰模型中的80 mm × 400 mm细长比例，
% 这里默认设置为80 mm × 80 mm正方形标定板。
obj_w = 80;           % 标定板横向宽度，mm
obj_h = 80;           % 标定板纵向高度，mm

fprintf('\n黑白标定板物理参数：\n');
fprintf('标定板宽度 obj_w = %.2f mm\n', obj_w);
fprintf('标定板高度 obj_h = %.2f mm\n', obj_h);
fprintf('标定板物距 d = %.2f mm\n', d);

%% ============================================================
%  5. 光场相机参数
%% ============================================================

v = 16;               % 主透镜到微透镜阵列距离，mm
N_line = 81;          % 主透镜采样数

D = 4;                % 主透镜直径，mm
F = 16;               % 主透镜焦距，mm

lens_d = 0.01;        % 微透镜直径，mm
sen_N = 20;           % 每个微透镜对应的传感器像素数

% 微透镜焦距
lens_f = lens_d * F / D;

% 微透镜阵列到传感器的相关距离
lens_v = lens_f + v;

% 微透镜数量
micr_N = fix(D / lens_d);

% 传感器总像素数
sen_N_total = micr_N * sen_N;

% 单个传感器像素尺寸
sen_d = lens_d / sen_N;

fprintf('\n光场相机参数：\n');
fprintf('主透镜直径 D = %.4f mm\n', D);
fprintf('主透镜焦距 F = %.4f mm\n', F);
fprintf('微透镜直径 lens_d = %.6f mm\n', lens_d);
fprintf('微透镜焦距 lens_f = %.6f mm\n', lens_f);
fprintf('微透镜数量 micr_N = %d\n', micr_N);
fprintf('每个微透镜后像素数 sen_N = %d\n', sen_N);
fprintf('传感器理论像素数 = %d × %d\n', ...
    sen_N_total, sen_N_total);

%% ============================================================
%  6. 将黑白标定板送入LF_sim
%
%  为兼容原来的reshape4to2_im，仍然分别保存RGB三个通道。
%% ============================================================

im_max = zeros(3, 1);

tic;

for k = 1:3

    fprintf('\n正在计算标定板第%d个通道的光场成像……\n', k);

    obje = double(obj(:, :, k));

    im = LF_sim( ...
        obje, ...
        d, ...
        obj_w, ...
        obj_h, ...
        D, ...
        N_line, ...
        F, ...
        v, ...
        lens_d, ...
        lens_f, ...
        sen_N);

    im_max(k) = max(im(:));

    output_mat_file = sprintf( ...
        './dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
        d, v, N_line, k);

    save(output_mat_file, 'im', '-v7.3');

    fprintf('第%d通道光场已保存：%s\n', ...
        k, output_mat_file);

    fprintf('第%d通道最大光场强度：%.6e\n', ...
        k, im_max(k));
end

simulation_time = toc;

fprintf('\n黑白标定板光场传输时间：%.3f s\n', ...
    simulation_time);

%% ============================================================
%  7. 加载RGB四维光场并组装为5D光场
%
%  LF_sim输出维度：
%  [空间x, 空间y, 角域u, 角域v]
%
%  组装后：
%  [空间x, 空间y, 角域u, 角域v, RGB]
%% ============================================================

LF_R_data = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_1.mat', ...
    d, v, N_line), 'im');

LF_G_data = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_2.mat', ...
    d, v, N_line), 'im');

LF_B_data = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_3.mat', ...
    d, v, N_line), 'im');

LF_R = LF_R_data.im;
LF_G = LF_G_data.im;
LF_B = LF_B_data.im;

fprintf('\nR通道四维光场尺寸：\n');
disp(size(LF_R));

fprintf('G通道四维光场尺寸：\n');
disp(size(LF_G));

fprintf('B通道四维光场尺寸：\n');
disp(size(LF_B));

if ~isequal(size(LF_R), size(LF_G), size(LF_B))
    error('RGB三个通道的四维光场尺寸不一致。');
end

% 组装5D光场
LF = cat(5, LF_R, LF_G, LF_B);

fprintf('\n黑白标定板5D光场尺寸：\n');
disp(size(LF));

fprintf('LF维数：%d\n', ndims(LF));

LF_filename = sprintf( ...
    './dataRGB/LF_5D_stepboard_d_%d_v_%d_Nline_%d.mat', ...
    d, v, N_line);

save(LF_filename, 'LF', '-v7.3');

fprintf('黑白标定板5D光场已保存：\n%s\n', ...
    LF_filename);

%% ============================================================
%  8. 生成二维传感器图像
%% ============================================================

im_RGB = reshape4to2_im( ...
    d, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure('Color', 'white');

imshow(im_RGB, []);

title('黑白阶跃标定板的微透镜阵列传感器图像');

sensor_filename = sprintf( ...
    './dataRGB/stepboard_sensor_d_%d_v_%d_Nline_%d_RGB.png', ...
    d, v, N_line);

imwrite(im_RGB, sensor_filename);

fprintf('\n传感器图像已保存：\n%s\n', ...
    sensor_filename);

%% ===================== 8. 原始数据求和重构 =====================

% 保留double精度，用于后续灰度分析和一阶差分
im_micr_raw = sum_4D_im(d, v, N_line, micr_N);
im_micr_raw = double(im_micr_raw);

fprintf('\n原始求和重构数据尺寸：\n');
disp(size(im_micr_raw));

% 如果函数返回RGB图像，则先转换为灰度图像
if ndims(im_micr_raw) == 3

    if size(im_micr_raw, 3) >= 3
        im_micr_gray = ...
            0.2989 .* im_micr_raw(:, :, 1) + ...
            0.5870 .* im_micr_raw(:, :, 2) + ...
            0.1140 .* im_micr_raw(:, :, 3);
    else
        im_micr_gray = im_micr_raw(:, :, 1);
    end

else
    im_micr_gray = im_micr_raw;
end

% 去除异常数值
im_micr_gray(~isfinite(im_micr_gray)) = 0;

% 全图归一化
im_min = min(im_micr_gray(:));
im_max = max(im_micr_gray(:));

if im_max > im_min
    im_micr_norm = ...
        (im_micr_gray - im_min) ./ (im_max - im_min);
else
    im_micr_norm = zeros(size(im_micr_gray));
end

% 生成用于显示和保存的uint8图像
im_micr_show = im2uint8(im_micr_norm);

figure('Color', 'white');
imshow(im_micr_show, []);
title('黑白阶跃标定板原始数据求和重构');

fprintf('\n已画出原始数据求和重构图像！\n');

