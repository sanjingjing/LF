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

%% ===================== 1. 分层火焰光场成像模型 =====================
Ny = 500;
Nz = 500;

[y, z] = meshgrid(linspace(-0.04, 0.04, Ny), ...
                  linspace(0, 0.4, Nz));

% 火焰温度模型参数
T_base = 1200;       % K
T_amp  = 600;        % K
m = 3;
n = 0.9;

R = 0.04;            % 火焰径向尺度 m
Z = 0.4;             % 火焰轴向尺度 m

% 光谱辐射参数
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K

% 衰减系数，按图中给出的预测衰减系数
ke = 10;             % m^-1

% 分层参数：x 方向从 -0.03 m 到 0.03 m，每层间隔 0.01 m
x_centers = -0.03 : 0.01 : 0.03;
n_layer = length(x_centers);

% 每层厚度范围为 x_center ± 0.005 m
layer_half_thickness = 0.005;     % m

% 火焰中心线距离主镜头 0.8 m
d_center = 0.8;                   % m

% 各层到主镜头距离，单位 mm
d_list = (d_center + x_centers) * 1000;



fprintf('各层到主镜头距离 d_list = \n');
disp(d_list);

%% ===================== 2. 体辐射积分与灰体衰减 =====================
% I_layers 保存每一层的出射面辐射图像
I_layers = zeros(Nz, Ny, n_layer);

% 为了模拟每层厚度内大量颗粒的体辐射，层内继续离散积分
sub_N = 21;
dx = 2 * layer_half_thickness / (sub_N - 1);   % m

% 灰体局部发射率 / 吸收率
% alpha = epsilon = 1 - exp(-ke * dx)
epsilon_dx = 1 - exp(-ke * dx);

for i = 1:n_layer

    x0 = x_centers(i);
    x_samples = linspace(x0 - layer_half_thickness, ...
                         x0 + layer_half_thickness, sub_N);

    Es_out = zeros(Nz, Ny);    % 当前层最终出射面辐射力

    for j = 1:sub_N

        x = x_samples(j);

        % 三维半径 r = sqrt(x^2 + y^2)
        r = sqrt(x.^2 + y.^2);

        % 图中给出的复杂温度场
        term = m * ((r.^2 ./ R.^2) + (z.^2 ./ Z.^2)) - n;
        T = T_base + T_amp .* exp(-(term.^2));

        % 黑体单色辐射
        Eb = c1 .* lambda.^(-5) ./ ...
            (exp(c2 ./ (lambda .* T)) - 1);

        % 当前体元到本层出射面的距离
        % 假设相机在 x 负方向，辐射向 x 负方向出射
        L_exit = x - (x0 - layer_half_thickness);

        % 灰体发射 + 衰减
        % epsilon = alpha = 1 - exp(-ke * dx)
        % 出射贡献：epsilon * Eb * exp(-ke * L_exit)
        Es_out = Es_out + epsilon_dx .* Eb .* exp(-ke .* L_exit);

    end

    % 这里 Es_out 是该层沿拍摄方向积分后的面辐射图像
    I_layers(:,:,i) = Es_out;
end

I_total = sum(I_layers, 3);

%% ===================== 3. 辐射强度转 RGB 火焰颜色 =====================


RGB_layers = zeros(Nz, Ny, 3, n_layer);

for i = 1:n_layer

    I_norm = I_layers(:,:,i);
    I_norm = I_norm - min(I_norm(:));

    if max(I_norm(:)) > 0
        I_norm = I_norm ./ max(I_norm(:));
    end

    % Gamma 调整，使火焰更明显
    I_norm = I_norm .^ 0.65;

    Rch = zeros(Nz, Ny);
    Gch = zeros(Nz, Ny);
    Bch = zeros(Nz, Ny);

    % 火焰 RGB 映射：黑 -> 红 -> 橙 -> 黄白
    Rch = min(1, 2.4 .* I_norm);
    Gch = min(1, 1.6 .* max(I_norm - 0.12, 0) ./ 0.88);
    Bch = min(1, 0.45 .* max(I_norm - 0.65, 0) ./ 0.35);

    RGB_layers(:,:,1,i) = Rch;
    RGB_layers(:,:,2,i) = Gch;
    RGB_layers(:,:,3,i) = Bch;
end

RGB_total = sum(RGB_layers, 4);

% 叠加后直接限幅到 [0,1]，不要再除以全局最大值
RGB_total = min(RGB_total, 1);

%% ===================== 4. 显示各层 RGB 火焰 =====================
figure('Color','white');

for i = 1:n_layer
    subplot(1, n_layer, i);

    RGB_show = RGB_layers(:,:,:,i);

    if max(RGB_show(:)) > 0
        RGB_show = RGB_show ./ max(RGB_show(:));
    end

    image(linspace(-0.04, 0.04, Ny), ...
          linspace(0, 0.4, Nz), ...
          RGB_show);

    set(gca, 'YDir', 'normal');
    axis image off;

    title(['第', num2str(i), '层 d=', num2str(d_list(i)), ' mm']);
end

%% ===================== 5. 显示总 RGB 火焰辐射图 =====================
figure('Color','white');

image(linspace(-0.04, 0.04, Ny), ...
      linspace(0, 0.4, Nz), ...
      RGB_total);

set(gca, 'YDir', 'normal');
axis image;
axis tight;
grid off;

xlabel('径向坐标 y (m)');
ylabel('轴向坐标 z (m)');
title('七层火焰叠加 RGB 辐射图');

%% ===================== 6. 归一化成 LF_sim 输入 =====================
obj_layers = zeros(Nz, Ny, 3, n_layer, 'uint8');

for i = 1:n_layer
    RGB_i = RGB_layers(:,:,:,i);

    if max(RGB_i(:)) > 0
        RGB_i = RGB_i ./ max(RGB_i(:));
    end

    obj_layers(:,:,:,i) = uint8(RGB_i * 255);
end

%% ===================== 7. 光场相机参数 =====================
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

obj_w = 0.08 * 1000;   % mm
obj_h = 0.4 * 1000;   % mm

%% ===================== 8. 分层传播并叠加 =====================
im_total_R = [];
im_total_G = [];
im_total_B = [];

tic

for k = 1:3

    fprintf('正在计算第 %d 个 RGB 通道...\n', k);

    im_sum = 0;

    for i = 1:n_layer

        fprintf('    正在传播第 %d 层火焰，距离 d = %.1f mm...\n', ...
            i, d_list(i));

        obje = double(obj_layers(:,:,k,i));

        im_i = LF_sim(obje, d_list(i), obj_w, obj_h, ...
            D, N_line, F, v, lens_d, lens_f, sen_N);

        im_sum = im_sum + im_i;

    end

    save(sprintf('./dataRGB/im_flame_layered_v_%d_Nline_%d_%d.mat', ...
        v, N_line, k), 'im_sum');

    if k == 1
        im_total_R = im_sum;
    elseif k == 2
        im_total_G = im_sum;
    else
        im_total_B = im_sum;
    end

end

t = toc;
fprintf('\n模拟分层 RGB 火焰光场传输时间 t = %.3f s\n', t);

%% ===================== 9. 组装 5D 光场 =====================
if ~isequal(size(im_total_R), size(im_total_G), size(im_total_B))
    error('三个通道尺寸不一致，无法组装 5D 光场');
end

LF = cat(5, im_total_R, im_total_G, im_total_B);

fprintf('\n5D 光场 LF 的维度为: ');
disp(size(LF));
fprintf('LF 的维数 ndims(LF) = %d\n', ndims(LF));

save(sprintf('./dataRGB/LF_5D_layered_flame_v_%d_Nline_%d.mat', ...
    v, N_line), 'LF', '-v7.3');

fprintf('5D 光场 LF 已保存完成！\n');

%% ===================== 10. 显示输入 LF_sim 的 RGB 火焰总图 =====================
obj_total = zeros(Nz, Ny, 3);

for i = 1:n_layer
    obj_total = obj_total + double(obj_layers(:,:,:,i));
end

if max(obj_total(:)) > 0
    obj_total = uint8(obj_total ./ max(obj_total(:)) * 255);
else
    obj_total = uint8(obj_total);
end

figure;
imshow(obj_total, []);
title('输入 LF\_sim 的分层 RGB 火焰总图');

%% ===================== 11. 显示微透镜阵列传感器 RGB 图像 =====================
im_max = [max(im_total_R(:)); ...
          max(im_total_G(:)); ...
          max(im_total_B(:))];

% 为了兼容 reshape4to2_im 的文件名读取格式
% 这里用中心层距离 d = 800 mm
d = 800;

im = im_total_R;
save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d, v, N_line, 1), 'im');

im = im_total_G;
save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d, v, N_line, 2), 'im');

im = im_total_B;
save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d, v, N_line, 3), 'im');

im_RGB = reshape4to2_im(d, v, N_line, micr_N, sen_N, im_max);

figure;
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像：分层 RGB 火焰');

imwrite(im_RGB, sprintf('./dataRGB/im_layered_flame_v_%d_Nline_%d_RGB.jpg', ...
    v, N_line));

fprintf('\n仿真完成：已考虑灰体发射、衰减和七层火焰深度分布。\n');


%% ===================== 12. 原始数据求和重构后的图像 =====================
im_micr = sum_4D_im(d, v, N_line, micr_N);

% 转成 double，便于显示增强
im_micr = double(im_micr);

% 不用最大值，改用 99% 分位数作为显示上限
max_im_micr = prctile(im_micr(:), 99);

if max_im_micr > 0
    im_micr = im_micr ./ max_im_micr;
end

% 饱和截断，避免整体被压暗
im_micr = min(im_micr, 1);

% gamma 增强，颜色更亮、更浓
im_micr = im_micr .^ 0.65;

% 转回 uint8
im_micr = uint8(im_micr * 255);

figure;
imshow(im_micr, []);
title('原始数据求和重构：分层 RGB 火焰');

fprintf('\n已画出原始数据重构图像!\n');
disp('The end!');