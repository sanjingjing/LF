clear all;
clc;
close all;

hsize = 37;
num_iter = 15;
image_size = 512;

% 是否显示时上下反转火焰
% true：显示为正像
% false：保持原始矩阵方向
flip_flame_show = true;

%% ===================== 1. 圆柱形火焰九层固定 Z_w 截面模型 =====================
t1 = 1200;        % 温度变化幅度 K
t2 = 900;         % 温度下限 K
m  = 3;           % 温度梯度参数
n  = 0.9;         % 高温层位置参数

RF = 0.02;        % 火焰半径 m
HF = 0.20;        % 火焰高度 m

% 九个固定 Z_w 截面位置
% 必须满足 abs(Z_layer) <= RF
Z_layer = [-0.009, -0.007, -0.005, -0.002, 0, ...
            0.002,  0.005,  0.007,  0.009];

layer_num = numel(Z_layer);

%% ===================== 2. 建立二维截面坐标 =====================
N = image_size;

X_range = linspace(-HF/2, HF/2, N);   % X_w 横向
Y_range = linspace(0, HF, N);         % Y_w 轴向高度

[Xw, Yw] = meshgrid(X_range, Y_range);

%% ===================== 3. 单色普朗克辐射模型参数 =====================
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K

%% ===================== 4. Beer 衰减参数 =====================
ke = 8;              % m^-1

% 假设相机位于 +Z_w 方向
Z_cam = RF;

%% ===================== 5. 生成九层火焰温度场与衰减辐射场 =====================
T_stack  = zeros(N, N, layer_num);
Eb_stack = zeros(N, N, layer_num);
Q_stack  = zeros(N, N, layer_num);
mask_stack = false(N, N, layer_num);
L_layer = zeros(layer_num, 1);

for i = 1:layer_num

    Zi = Z_layer(i);

    % 当前 Z_w 截面的圆柱边界
    r2 = Xw.^2 + Zi.^2;
    mask = (r2 <= RF^2) & (Yw >= 0) & (Yw <= HF);

    % 温度场模型
    term = m .* (Xw.^2 ./ RF.^2 + Zi.^2 ./ RF.^2 + Yw.^2 ./ HF.^2) - n;
    T = t1 .* exp(-(term.^2)) + t2;
    T(~mask) = t2;

    % 单色普朗克辐射
    Eb = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T)) - 1);
    Eb(~mask) = 0;

    % Beer 衰减
    L = Z_cam - Zi;
    Q = Eb .* exp(-ke .* L);

    T_stack(:,:,i) = T;
    Eb_stack(:,:,i) = Eb;
    Q_stack(:,:,i) = Q;
    mask_stack(:,:,i) = mask;
    L_layer(i) = L;
end

%% ===================== 6. 显示九层真实温度场 =====================
figure('Color','white');
for i = 1:layer_num
    subplot(3,3,i);
    imshow(orientForShow(T_stack(:,:,i), flip_flame_show), []);
    colormap hot;
    colorbar;
    title(['T_', num2str(i), ', Z=', num2str(Z_layer(i)), ' m']);
end
sgtitle('九层真实温度场，仅用于验证火焰模型');

%% ===================== 7. 显示九层衰减后辐射图像 =====================
figure('Color','white');
for i = 1:layer_num
    subplot(3,3,i);
    imshow(orientForShow(Q_stack(:,:,i), flip_flame_show), []);
    title(['Q_', num2str(i)]);
end
sgtitle('九层衰减后辐射图像');

I_direct = sum(Q_stack, 3);

figure('Color','white');
imshow(orientForShow(I_direct, flip_flame_show), []);
title('九层火焰辐射直接叠加图像');

%% ===================== 8. 构造 9×9 点扩散函数矩阵 =====================
% Op{k,i} 表示第 i 层火焰对第 k 张重聚焦图像的贡献
% k = i 时表示聚焦层，sigma 小；
% k ~= i 时表示离焦层，sigma 大。
%
% 前 3×3 保留你原始三层模型中的点扩散参数：
% [2.29 6.81 10.73;
%  6.36 2.17 6.59;
%  9.05 6.19 2.61]
%
% 其余项按照“离焦距离越远，sigma 越大”的原则补充。

PSF_sigma = [
    2.29,  6.81, 10.73, 12.20, 13.40, 14.20, 15.00, 15.60, 16.00;
    6.36,  2.17,  6.59, 10.20, 12.00, 13.20, 14.10, 15.00, 15.60;
    9.05,  6.19,  2.61,  6.70, 10.30, 12.10, 13.30, 14.20, 15.00;
   11.20,  9.10,  6.36,  2.29,  6.81, 10.73, 12.20, 13.40, 14.20;
   12.40, 10.20,  6.36,  6.36,  2.17,  6.59, 10.20, 12.00, 13.20;
   13.20, 11.30,  9.05,  6.19,  6.19,  2.61,  6.70, 10.30, 12.10;
   14.00, 12.40, 10.60,  9.05,  6.36,  6.81,  2.29,  6.81, 10.73;
   14.70, 13.20, 11.60, 10.20,  9.10,  6.36,  6.36,  2.17,  6.59;
   15.40, 14.00, 12.70, 11.40, 10.00,  9.05,  9.05,  6.19,  2.61
];

Op = cell(layer_num, layer_num);
Op_adj = cell(layer_num, layer_num);

for k = 1:layer_num
    for i = 1:layer_num

        sigma = PSF_sigma(k,i);

        [Op{k,i}, Op_adj{k,i}, ~] = ...
            createGaussianBlurringOperator([N N], hsize, sigma);
    end
end

disp('9×9 点扩散函数 sigma 矩阵为：');
disp(PSF_sigma);

figure('Color','white');
imagesc(PSF_sigma);
axis image;
colorbar;
title('9×9 点扩散函数 sigma 矩阵');
xlabel('火焰层 i');
ylabel('重聚焦图像 k');

%% ===================== 9. Forward imaging model：生成九张重聚焦叠加图像 =====================
% I_focus(:,:,k) 为重聚焦到第 k 层时的图像
% 每张图像都是 9 层火焰经过不同 PSF 后的叠加

I_focus = zeros(N, N, layer_num);

for k = 1:layer_num
    temp = zeros(N, N);

    for i = 1:layer_num
        temp = temp + Op{k,i}(Q_stack(:,:,i));
    end

    I_focus(:,:,k) = temp;
end

figure('Color','white');
for k = 1:layer_num
    subplot(3,3,k);
    imshow(orientForShow(I_focus(:,:,k), flip_flame_show), []);
    title(['重聚焦叠加图像 I_', num2str(k)]);
end
sgtitle('九张重聚焦叠加图像');

%% ===================== 10. Stack measured images b =====================
% image = [I1; I2; ...; I9]

image = zeros(layer_num * N, N);

for k = 1:layer_num
    row_start = (k-1)*N + 1;
    row_end   = k*N;

    image(row_start:row_end, :) = I_focus(:,:,k);
end

image = max(real(image), 0);

% 仅用于显示，把每一块单独上下翻转
image_show = zeros(size(image));

for k = 1:layer_num
    row_start = (k-1)*N + 1;
    row_end   = k*N;

    image_show(row_start:row_end, :) = ...
        orientForShow(image(row_start:row_end, :), flip_flame_show);
end

figure('Color','white');
imshow(image_show, []);
title('观测叠加图像 b = [I_1; I_2; ...; I_9]');

%% ===================== 11. 定义九层块前向算子 A =====================
A = @(x) applyA_9layer(x, Op, N, layer_num);

%% ===================== 12. 定义九层伴随算子 A_adj =====================
A_adj = @(y) applyAT_9layer(y, Op_adj, N, layer_num);

%% ===================== 13. Richardson-Lucy deconvolution =====================
% 初始估计必须为正
init = ones(size(image)) .* max(mean(image(:)), eps);

% RL 重建，得到九层衰减后的辐射强度估计值 Q_i
field_rl = deconvRL_operator(A, A_adj, image, num_iter, init);

field_rl = real(field_rl);
field_rl(field_rl < 0) = 0;

%% ===================== 14. 拆分九层 RL 恢复结果 =====================
recQ_stack = zeros(N, N, layer_num);

for i = 1:layer_num
    row_start = (i-1)*N + 1;
    row_end   = i*N;

    recQ_stack(:,:,i) = field_rl(row_start:row_end, :);
end

%% ===================== 15. 显示 RL 重建整体结果：拆成 3×3 =====================
figure('Color','white');

for i = 1:layer_num
    row_start = (i-1)*N + 1;
    row_end   = i*N;

    rec_i = field_rl(row_start:row_end, :);

    subplot(3,3,i);
    imshow(orientForShow(rec_i, flip_flame_show), []);
    title(['RL 恢复 Q_', num2str(i)]);
end

sgtitle('RL 去卷积整体结果 3×3 显示');

%% ===================== 16. 显示九层 RL 重建结果：统一灰度范围 =====================

recQ_show_max = max(recQ_stack(:));
recQ_show_min = min(recQ_stack(:));

figure('Color','white');

for i = 1:layer_num
    subplot(3,3,i);

    temp = recQ_stack(:,:,i);

    if flip_flame_show
        temp = flipud(temp);
    end

    imshow(temp, [recQ_show_min, recQ_show_max]);
    title(['RL 恢复 Q_', num2str(i)]);
end

sgtitle('九层 RL 重建后的衰减辐射图像：统一灰度范围显示');

%% ===================== 本脚本用到的局部函数 =====================

function y = applyA_9layer(x, Op, N, layer_num)
% 前向算子：
% 输入 x 为 [Q1; Q2; ...; Q9]
% 输出 y 为 [I1; I2; ...; I9]

    y = zeros(layer_num * N, N);

    for k = 1:layer_num
        temp = zeros(N, N);

        for i = 1:layer_num
            row_start = (i-1)*N + 1;
            row_end   = i*N;

            Qi = x(row_start:row_end, :);

            temp = temp + Op{k,i}(Qi);
        end

        out_start = (k-1)*N + 1;
        out_end   = k*N;

        y(out_start:out_end, :) = temp;
    end
end

function x = applyAT_9layer(y, Op_adj, N, layer_num)
% 伴随算子：
% 输入 y 为 [I1; I2; ...; I9]
% 输出 x 为 [Q1; Q2; ...; Q9]

    x = zeros(layer_num * N, N);

    for i = 1:layer_num
        temp = zeros(N, N);

        for k = 1:layer_num
            row_start = (k-1)*N + 1;
            row_end   = k*N;

            Ik = y(row_start:row_end, :);

            temp = temp + Op_adj{k,i}(Ik);
        end

        out_start = (i-1)*N + 1;
        out_end   = i*N;

        x(out_start:out_end, :) = temp;
    end
end

function out = orientForShow(A, flip_flag)
% 只用于显示方向修正，不改变计算矩阵本身

    out = mat2gray(A);

    if flip_flag
        out = flipud(out);
    end
end