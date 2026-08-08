clear; clc; close all;

addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\频域重聚焦精简'));

%% ===================== 1. 火焰物理区域 =====================
t1 = 1200;
t2 = 900;
m  = 3;
n  = 0.9;
R0 = 0.04;      % 火焰特征径向尺度
Z0 = 0.4;       % 火焰特征轴向尺度
x0 = 0;

Ny = 500;
Nz = 500;

% ===== 关键：扩大计算域，而不是人为加黑边 =====
r_min = -0.04;
r_max =  0.04;
z_min =  0.00;
z_max =  0.40;

[y,z] = meshgrid(linspace(r_min,r_max,Ny), linspace(z_min,z_max,Nz));

%% ===================== 2. 温度场 =====================
term = m * (x0.^2/R0.^2 + y.^2/R0.^2 + z.^2/Z0.^2) - n;
T = t1 * exp(-(term.^2)) + t2;

%% ===================== 3. 单色辐射模型 =====================
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K
Eb = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T)) - 1);

hot_mask = Eb > 0.02 * max(Eb(:));

%% ===================== 4. Monte Carlo 辐射传输 =====================
mc_params.ke = 30;
mc_params.omega = 0.0;
mc_params.N_rays = 80;
mc_params.max_scatter = 6;
mc_params.use_cone_emission = true;
mc_params.alpha_deg = 0.6;
mc_params.use_anisotropic = false;
mc_params.A1 = 0.0;

mc_params.ymin = min(y(:));
mc_params.ymax = max(y(:));
mc_params.zmin = min(z(:));
mc_params.zmax = max(z(:));

fprintf('开始 Monte Carlo 辐射传输...\n');
tic;
I_mc = mc_radiative_transfer_2d(Eb, hot_mask, y, z, mc_params);
toc;

%% ===================== 5. 显示 MC 火焰结果 =====================
figure('Color','w');
imagesc(linspace(r_min,r_max,Ny), linspace(z_min,z_max,Nz), I_mc);
set(gca,'YDir','normal');
colormap(gray(256));
colorbar;
axis image;
pbaspect([1 5 1]);
xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');
title('Monte Carlo 火焰辐射传输结果');

%% ===================== 6. 直接用 MC 结果构造 LF_sim 输入 =====================
% 不再创建黑色背景，不再人为贴图
I_norm = I_mc - min(I_mc(:));
if max(I_norm(:)) > 0
    I_norm = I_norm / max(I_norm(:));
end

% 显示压缩，仅用于数值映射到 8bit，不改变物理区域
gamma_val = 0.75;
I_show = I_norm .^ gamma_val;

obj_gray = uint8(I_show * 255);
obj = repmat(obj_gray, [1,1,3]);

figure('Color','w');
imshow(obj, []);
title('直接由 MC 火焰辐射场构造的 LF\_sim 输入');

%% ===================== 7. 光场相机参数 =====================
v = 16;
N_line = 81;

D = 4;
F = 16;
lens_d = 0.08;
sen_N = 20;

lens_f = lens_d * F / D;
micr_N = ceil(D / lens_d);

im_max = zeros(3,1);
d = 600;

%% ===================== 8. LF_sim 成像 =====================
tic;
for k = 1:3
    fprintf('正在计算第 %d 个通道...\n', k);
    obje = double(obj(:,:,k));

    % ===== 关键：物面尺寸直接对应真实物理区域 =====
    obj_w = 0.15 * 1000;   % mm
    obj_h = 0.45 * 1000;   % mm

    im = LF_sim(obje, d, obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im_max(k) = max(im(:));

    save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, k), 'im');
end
toc;

%% ===================== 9. 显示微透镜传感器图像 =====================
im_RGB = reshape4to2_im(d, v, N_line, micr_N, sen_N, im_max);

figure('Color','w');
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像');
axis image;
pbaspect([1 5 1]);

imwrite(im_RGB, sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_RGB.jpg', d, v, N_line));

disp('The end!');
