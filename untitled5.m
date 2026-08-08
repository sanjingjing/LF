clear all;
clc;
close all;

%% ===================== 参数设置 =====================
hsize = 37;
num_iter = 10;
image_size = 512;

%% ===================== 1. 火焰辐射模型 =====================
t1 = 1200;      % 温度增量 K
t2 = 900;       % 环境温度 K
m  = 10;        % 形状参数
n  = 0.4;
R  = 0.025;     % m
Z  = 0.075;     % m

x1 = 0;
x2 = 0.005;
x3 = 0.008;

Ny = 500;
Nz = 500;

[y,z] = meshgrid( ...
    linspace(-0.02, 0.02, Ny), ...
    linspace(0, 0.075, Nz));

%% ===================== 2. 温度场 =====================
term1 = m * (x1.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term2 = m * (x2.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;
term3 = m * (x3.^2/R.^2 + y.^2/R.^2 + z.^2/Z.^2) - n;

T1 = t1 * exp(-(term1.^2)) + t2;
T2 = t1 * exp(-(term2.^2)) + t2;
T3 = t1 * exp(-(term3.^2)) + t2;

figure;
subplot(1,3,1); imagesc(y(1,:), z(:,1), T1); axis image; colormap hot; colorbar; title('温度场 T1');
subplot(1,3,2); imagesc(y(1,:), z(:,1), T2); axis image; colormap hot; colorbar; title('温度场 T2');
subplot(1,3,3); imagesc(y(1,:), z(:,1), T3); axis image; colormap hot; colorbar; title('温度场 T3');

%% ===================== 3. RGB 单色辐射模型 =====================
% 单位：um
lambdaR = 0.65;
lambdaG = 0.55;
lambdaB = 0.45;

c1 = 3.7418e-16;
c2 = 1.4388e4;

Planck = @(T, lambda) c1 .* lambda.^(-5) ./ ...
    (exp(c2 ./ (lambda .* T)) - 1);

Eb1_R = Planck(T1, lambdaR);
Eb1_G = Planck(T1, lambdaG);
Eb1_B = Planck(T1, lambdaB);

Eb2_R = Planck(T2, lambdaR);
Eb2_G = Planck(T2, lambdaG);
Eb2_B = Planck(T2, lambdaB);

Eb3_R = Planck(T3, lambdaR);
Eb3_G = Planck(T3, lambdaG);
Eb3_B = Planck(T3, lambdaB);

%% ===================== 4. Beer 衰减 =====================
ke = 8;               % m^-1
x_max = 0.025;

L1 = x_max - x1;
L2 = x_max - x2;
L3 = x_max - x3;

Q1_R = Eb1_R .* exp(-ke .* L1);
Q1_G = Eb1_G .* exp(-ke .* L1);
Q1_B = Eb1_B .* exp(-ke .* L1);

Q2_R = Eb2_R .* exp(-ke .* L2);
Q2_G = Eb2_G .* exp(-ke .* L2);
Q2_B = Eb2_B .* exp(-ke .* L2);

Q3_R = Eb3_R .* exp(-ke .* L3);
Q3_G = Eb3_G .* exp(-ke .* L3);
Q3_B = Eb3_B .* exp(-ke .* L3);

%% ===================== 5. resize + 翻转 =====================
Q1_R = flipud(imresize(im2double(Q1_R), [image_size image_size]));
Q1_G = flipud(imresize(im2double(Q1_G), [image_size image_size]));
Q1_B = flipud(imresize(im2double(Q1_B), [image_size image_size]));

Q2_R = flipud(imresize(im2double(Q2_R), [image_size image_size]));
Q2_G = flipud(imresize(im2double(Q2_G), [image_size image_size]));
Q2_B = flipud(imresize(im2double(Q2_B), [image_size image_size]));

Q3_R = flipud(imresize(im2double(Q3_R), [image_size image_size]));
Q3_G = flipud(imresize(im2double(Q3_G), [image_size image_size]));
Q3_B = flipud(imresize(im2double(Q3_B), [image_size image_size]));

%% ===================== 6. 组合成 RGB 三层真实原始亮度 =====================
Q1_rgb = zeros(image_size, image_size, 3);
Q2_rgb = zeros(image_size, image_size, 3);
Q3_rgb = zeros(image_size, image_size, 3);

Q1_rgb(:,:,1) = Q1_R;
Q1_rgb(:,:,2) = Q1_G;
Q1_rgb(:,:,3) = Q1_B;

Q2_rgb(:,:,1) = Q2_R;
Q2_rgb(:,:,2) = Q2_G;
Q2_rgb(:,:,3) = Q2_B;

Q3_rgb(:,:,1) = Q3_R;
Q3_rgb(:,:,2) = Q3_G;
Q3_rgb(:,:,3) = Q3_B;

%% ===================== 7. RGB 强度归一化，保持层间相对关系 =====================
all_rgb = cat(4, Q1_rgb, Q2_rgb, Q3_rgb);

for c = 1:3
    max_c = max(all_rgb(:,:,c,:), [], 'all');
    Q1_rgb(:,:,c) = Q1_rgb(:,:,c) ./ max_c;
    Q2_rgb(:,:,c) = Q2_rgb(:,:,c) ./ max_c;
    Q3_rgb(:,:,c) = Q3_rgb(:,:,c) ./ max_c;
end

figure;
subplot(1,3,1); imshow(Q1_rgb); title('RGB 原始层 Q1');
subplot(1,3,2); imshow(Q2_rgb); title('RGB 原始层 Q2');
subplot(1,3,3); imshow(Q3_rgb); title('RGB 原始层 Q3');

%% ===================== 8. 原始三层叠加 RGB 图像 =====================
I_rgb = Q1_rgb + Q2_rgb + Q3_rgb;
I_rgb = I_rgb ./ max(I_rgb(:));

figure;
imshow(I_rgb);
title('三层火焰直接叠加 RGB 图像');

%% ===================== 9. 创建 Gaussian PSF 模糊算子 =====================
[Op11, Op11_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 2.29);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 6.81);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 10.73);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 6.36);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 2.17);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 6.59);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 9.05);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 6.19);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator([image_size image_size], hsize, 2.61);

%% ===================== 10. RGB 前向重聚焦成像 =====================
I1_rgb = zeros(image_size, image_size, 3);
I2_rgb = zeros(image_size, image_size, 3);
I3_rgb = zeros(image_size, image_size, 3);

for c = 1:3
    Q1 = Q1_rgb(:,:,c);
    Q2 = Q2_rgb(:,:,c);
    Q3 = Q3_rgb(:,:,c);

    I1_rgb(:,:,c) = Op11(Q1) + Op12(Q2) + Op13(Q3);
    I2_rgb(:,:,c) = Op21(Q1) + Op22(Q2) + Op23(Q3);
    I3_rgb(:,:,c) = Op31(Q1) + Op32(Q2) + Op33(Q3);
end

% 仅用于显示
show_I1 = I1_rgb ./ max(I1_rgb(:));
show_I2 = I2_rgb ./ max(I2_rgb(:));
show_I3 = I3_rgb ./ max(I3_rgb(:));

figure;
subplot(1,3,1); imshow(show_I1); title('RGB 重聚焦图像 1');
subplot(1,3,2); imshow(show_I2); title('RGB 重聚焦图像 2');
subplot(1,3,3); imshow(show_I3); title('RGB 重聚焦图像 3');

%% ===================== 11. 定义块矩阵前向算子 A =====================
A = @(x) [
    Op11(x(1:image_size, :)) + Op12(x(image_size+1:2*image_size, :)) + Op13(x(2*image_size+1:3*image_size, :));
    Op21(x(1:image_size, :)) + Op22(x(image_size+1:2*image_size, :)) + Op23(x(2*image_size+1:3*image_size, :));
    Op31(x(1:image_size, :)) + Op32(x(image_size+1:2*image_size, :)) + Op33(x(2*image_size+1:3*image_size, :))
];

%% ===================== 12. 定义伴随算子 A_adj =====================
A_adj = @(y) [
    Op11_adj(y(1:image_size, :)) + Op21_adj(y(image_size+1:2*image_size, :)) + Op31_adj(y(2*image_size+1:3*image_size, :));
    Op12_adj(y(1:image_size, :)) + Op22_adj(y(image_size+1:2*image_size, :)) + Op32_adj(y(2*image_size+1:3*image_size, :));
    Op13_adj(y(1:image_size, :)) + Op23_adj(y(image_size+1:2*image_size, :)) + Op33_adj(y(2*image_size+1:3*image_size, :))
];

%% ===================== 13. RGB 三通道分别 RL 去卷积 =====================
field_rl_rgb = zeros(3*image_size, image_size, 3);

for c = 1:3

    fprintf('正在进行第 %d 个 RGB 通道的 RL 去卷积...\n', c);

    image_c = zeros(3*image_size, image_size);

    image_c(1:image_size, :) = I1_rgb(:,:,c);
    image_c(image_size+1:2*image_size, :) = I2_rgb(:,:,c);
    image_c(2*image_size+1:3*image_size, :) = I3_rgb(:,:,c);

    image_c = max(real(image_c), 0);

    init_c = ones(size(image_c)) .* mean(image_c(:));

    field_c = deconvRL_operator(A, A_adj, image_c, num_iter, init_c);

    field_c = real(field_c);
    field_c(field_c < 0) = 0;

    field_rl_rgb(:,:,c) = field_c;
end

%% ===================== 14. 拆分恢复层 =====================
rec1_rgb = field_rl_rgb(1:image_size, :, :);
rec2_rgb = field_rl_rgb(image_size+1:2*image_size, :, :);
rec3_rgb = field_rl_rgb(2*image_size+1:3*image_size, :, :);

%% ===================== 15. 恢复结果统一归一化显示 =====================
rec_all = cat(4, rec1_rgb, rec2_rgb, rec3_rgb);

for c = 1:3
    max_c = max(rec_all(:,:,c,:), [], 'all');
    if max_c > 0
        rec1_rgb(:,:,c) = rec1_rgb(:,:,c) ./ max_c;
        rec2_rgb(:,:,c) = rec2_rgb(:,:,c) ./ max_c;
        rec3_rgb(:,:,c) = rec3_rgb(:,:,c) ./ max_c;
    end
end

figure;
subplot(1,3,1); imshow(rec1_rgb); title('RGB 恢复层 x1');
subplot(1,3,2); imshow(rec2_rgb); title('RGB 恢复层 x2');
subplot(1,3,3); imshow(rec3_rgb); title('RGB 恢复层 x3');

%% ===================== 16. 原始层与恢复层对比 =====================
figure;

subplot(2,3,1); imshow(Q1_rgb); title('真实层 Q1');
subplot(2,3,2); imshow(Q2_rgb); title('真实层 Q2');
subplot(2,3,3); imshow(Q3_rgb); title('真实层 Q3');

subplot(2,3,4); imshow(rec1_rgb); title('恢复层 x1');
subplot(2,3,5); imshow(rec2_rgb); title('恢复层 x2');
subplot(2,3,6); imshow(rec3_rgb); title('恢复层 x3');

disp('RGB 多通道 RL 去卷积完成。');

