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
[y, z] = meshgrid(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz));

% 温度参数
T_default = 0;
T_hot = 1373;

% 像素坐标网格
[Xp, Zp] = meshgrid(1:Ny, 1:Nz);

% 四个图形在图像中的中心位置（四个角）
cx1 = 150;  cz1 = 120;   % 左上：正方形
cx2 = 350;  cz2 = 120;   % 右上：圆形
cx3 = 150;  cz3 = 360;   % 左下：三角形
cx4 = 350;  cz4 = 360;   % 右下：正五边形

% 四层到主镜头距离（mm）
% 每层间隔 1 mm
d_list = [400,600, 800, 1000];
% 你也可以反过来用 [600, 601, 602, 603]
% 只要保持层间距离是 1 mm 即可

%% ===================== 2. 构造四层独立温度场 =====================
%% ===================== 2. 构造四层独立棋盘格温度场 =====================
T1 = ones(Nz, Ny) * T_default;
T2 = ones(Nz, Ny) * T_default;
T3 = ones(Nz, Ny) * T_default;
T4 = ones(Nz, Ny) * T_default;

% 棋盘格参数
checker_size = 100;   % 棋盘总边长，像素
num_square = 4;       % 4 × 4 棋盘
cell_size = checker_size / num_square;

% 生成棋盘格 mask 的函数
make_checker = @(cx,cz) ...
    checkerboard_mask(Xp, Zp, cx, cz, checker_size, num_square);

% 四个不同深度的棋盘格
mask1 = make_checker(cx1, cz1);   % 第1层：左上
mask2 = make_checker(cx2, cz2);   % 第2层：右上
mask3 = make_checker(cx3, cz3);   % 第3层：左下
mask4 = make_checker(cx4, cz4);   % 第4层：右下

% 黑白棋盘：白格为高温，黑格为背景
T1(mask1) = T_hot;
T2(mask2) = T_hot;
T3(mask3) = T_hot;
T4(mask4) = T_hot;

%% ===================== 3. 单色辐射模型 =====================
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K

Eb1 = zeros(size(T1));
Eb2 = zeros(size(T2));
Eb3 = zeros(size(T3));
Eb4 = zeros(size(T4));

idx1 = T1 > 0;
idx2 = T2 > 0;
idx3 = T3 > 0;
idx4 = T4 > 0;

Eb1(idx1) = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T1(idx1))) - 1);
Eb2(idx2) = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T2(idx2))) - 1);
Eb3(idx3) = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T3(idx3))) - 1);
Eb4(idx4) = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T4(idx4))) - 1);

%% ===================== 4. Beer 衰减 =====================
ke = 10;
y_max = max(y(:));
L = y_max - y;

I1 = Eb1 .* exp(-ke .* L);
I2 = Eb2 .* exp(-ke .* L);
I3 = Eb3 .* exp(-ke .* L);
I4 = Eb4 .* exp(-ke .* L);

% 总辐射图（仅用于显示）
I_total = I1 + I2 + I3 + I4;

%% ===================== 5. 显示各层棋盘格 =====================
shape1 = zeros(Nz, Ny); shape1(mask1) = 1;
shape2 = zeros(Nz, Ny); shape2(mask2) = 1;
shape3 = zeros(Nz, Ny); shape3(mask3) = 1;
shape4 = zeros(Nz, Ny); shape4(mask4) = 1;

figure('Color','white');
subplot(2,2,1); imshow(shape1,[]); title(['第1层 棋盘格, d = ', num2str(d_list(1)), ' mm']);
subplot(2,2,2); imshow(shape2,[]); title(['第2层 棋盘格, d = ', num2str(d_list(2)), ' mm']);
subplot(2,2,3); imshow(shape3,[]); title(['第3层 棋盘格, d = ', num2str(d_list(3)), ' mm']);
subplot(2,2,4); imshow(shape4,[]); title(['第4层 棋盘格, d = ', num2str(d_list(4)), ' mm']);

%% ===================== 6. 显示总辐射图 =====================
figure('Color','white');
imagesc(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz), I_total);
set(gca,'YDir','normal');
colormap(gray(256));
colorbar;
axis tight;
axis image;
grid off;
xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');
title('四层不同深度图形总辐射图');

%% ===================== 7. 归一化成 LF_sim 输入 =====================
I1_norm = I1 - min(I1(:)); if max(I1_norm(:))>0, I1_norm = I1_norm/max(I1_norm(:)); end
I2_norm = I2 - min(I2(:)); if max(I2_norm(:))>0, I2_norm = I2_norm/max(I2_norm(:)); end
I3_norm = I3 - min(I3(:)); if max(I3_norm(:))>0, I3_norm = I3_norm/max(I3_norm(:)); end
I4_norm = I4 - min(I4(:)); if max(I4_norm(:))>0, I4_norm = I4_norm/max(I4_norm(:)); end

I1_norm = I1_norm .^ 0.7;
I2_norm = I2_norm .^ 0.7;
I3_norm = I3_norm .^ 0.7;
I4_norm = I4_norm .^ 0.7;

obj1_gray = uint8(I1_norm * 255);
obj2_gray = uint8(I2_norm * 255);
obj3_gray = uint8(I3_norm * 255);
obj4_gray = uint8(I4_norm * 255);

obj1 = repmat(obj1_gray, [1,1,3]);
obj2 = repmat(obj2_gray, [1,1,3]);
obj3 = repmat(obj3_gray, [1,1,3]);
obj4 = repmat(obj4_gray, [1,1,3]);

%% ===================== 8. 光场相机参数 =====================

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
%% ===================== 本地函数：生成黑白棋盘格 =====================
function mask = checkerboard_mask(Xp, Zp, cx, cz, checker_size, num_square)

    cell_size = checker_size / num_square;

    x0 = cx - checker_size/2;
    z0 = cz - checker_size/2;

    x_rel = Xp - x0;
    z_rel = Zp - z0;

    inside = x_rel >= 0 & x_rel < checker_size & ...
             z_rel >= 0 & z_rel < checker_size;

    ix = floor(x_rel / cell_size);
    iz = floor(z_rel / cell_size);

    % 黑白交替，白格为 1
    pattern = mod(ix + iz, 2) == 0;

    mask = inside & pattern;
end