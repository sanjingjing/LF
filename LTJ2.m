clear all;
clc;

%% ===================== 添加路径 =====================
addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\火焰频域重聚焦'));

%% ===================== 1. 火焰辐射模型（四角分布 + 空心边框图形） =====================
Ny = 500;
Nz = 500;
[y, z] = meshgrid(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz));

% 初始化温度场
T_default = 0;
T_hot = 1373;
T = ones(Nz, Ny) * T_default;

% 像素坐标网格
[Xp, Zp] = meshgrid(1:Ny, 1:Nz);

% ---------------- 四个图形中心位置（四个角） ----------------
% 左上：正方形
cx1 = 150;  cz1 = 120;
% 右上：圆形
cx2 = 350;  cz2 = 120;
% 左下：三角形
cx3 = 150;  cz3 = 360;
% 右下：正五边形
cx4 = 350;  cz4 = 360;

%% ---------------- 第1层：正方形边框（左上） ----------------
side_outer = 90;   % 外边长
side_inner = 58;   % 内边长

mask1_outer = abs(Xp - cx1) <= side_outer/2 & abs(Zp - cz1) <= side_outer/2;
mask1_inner = abs(Xp - cx1) <= side_inner/2 & abs(Zp - cz1) <= side_inner/2;
mask1 = mask1_outer & ~mask1_inner;

%% ---------------- 第2层：圆环（右上） ----------------
r_outer = 48;
r_inner = 30;

dist2 = (Xp - cx2).^2 + (Zp - cz2).^2;
mask2 = (dist2 <= r_outer^2) & (dist2 >= r_inner^2);

%% ---------------- 第3层：三角形边框（左下，已放大） ----------------
tri_h_outer = 118;   % 外三角形高度
tri_w_outer = 126;   % 外三角形底边
tri_h_inner = 78;    % 内三角形高度
tri_w_inner = 82;    % 内三角形底边

% 外三角形：顶点朝上
x_tri_outer = [cx3, cx3 - tri_w_outer/2, cx3 + tri_w_outer/2];
z_tri_outer = [cz3 - tri_h_outer/2, cz3 + tri_h_outer/2, cz3 + tri_h_outer/2];

% 内三角形：同心缩小
x_tri_inner = [cx3, cx3 - tri_w_inner/2, cx3 + tri_w_inner/2];
z_tri_inner = [cz3 - tri_h_inner/2, cz3 + tri_h_inner/2, cz3 + tri_h_inner/2];

mask3_outer = inpolygon(Xp, Zp, x_tri_outer, z_tri_outer);
mask3_inner = inpolygon(Xp, Zp, x_tri_inner, z_tri_inner);
mask3 = mask3_outer & ~mask3_inner;

%% ---------------- 第4层：正五边形边框（右下） ----------------
R5_outer = 52;
R5_inner = 34;
theta0 = -pi/2;   % 让一个顶点朝上
theta = theta0 + (0:4) * 2*pi/5;

x_pent_outer = cx4 + R5_outer * cos(theta);
z_pent_outer = cz4 + R5_outer * sin(theta);

x_pent_inner = cx4 + R5_inner * cos(theta);
z_pent_inner = cz4 + R5_inner * sin(theta);

mask4_outer = inpolygon(Xp, Zp, x_pent_outer, z_pent_outer);
mask4_inner = inpolygon(Xp, Zp, x_pent_inner, z_pent_inner);
mask4 = mask4_outer & ~mask4_inner;

%% ---------------- 合成温度场 ----------------
T(mask1) = T_hot;
T(mask2) = T_hot;
T(mask3) = T_hot;
T(mask4) = T_hot;

%% ===================== 2. 单色辐射模型 =====================
lambda = 0.610;      % um
c1 = 3.7418e-16;
c2 = 1.4388e4;       % um*K

Eb = zeros(size(T));
validMask = T > 0;
Eb(validMask) = c1 .* lambda.^(-5) ./ (exp(c2 ./ (lambda .* T(validMask))) - 1);

%% ===================== 3. Beer 衰减 =====================
ke = 10;
y_max = max(y(:));
L = y_max - y;
I = Eb .* exp(-ke .* L);

%% ===================== 4. 显示火焰辐射图 =====================
figure('Color', 'white');
imagesc(linspace(-0.04, 0.04, Ny), linspace(0, 0.4, Nz), I);
set(gca, 'YDir', 'normal');
colormap(gray(256));

cb = colorbar;
cb.Label.String = '单色辐射强度';

axis tight;
axis image;
grid off;
xlabel('径向坐标 r (m)');
ylabel('轴向坐标 z (m)');
title(['四角分布几何图形辐射图, k_e = ', num2str(ke), ' m^{-1}']);

%% ===================== 5. 显示二值图形模型（便于检查） =====================
shape_img = zeros(Nz, Ny);
shape_img(mask1) = 1;
shape_img(mask2) = 1;
shape_img(mask3) = 1;
shape_img(mask4) = 1;

figure('Color', 'white');
imshow(shape_img, []);
title('四角分布空心图形：左上正方形，右上圆形，左下三角形，右下正五边形');

%% ===================== 6. 转为 LF_sim 输入 =====================
I_norm = I - min(I(:));
if max(I_norm(:)) > 0
    I_norm = I_norm / max(I_norm(:));
end

I_norm = I_norm .^ 0.7;

obj_gray = uint8(I_norm * 255);
obj = repmat(obj_gray, [1, 1, 3]);

d = 600;   % mm，发光体到主透镜距离

%% ===================== 7. 光场相机参数 =====================
v = 16;
N_line = 81;
D = 4;
F = 16;
lens_d = 0.04;
sen_N = 20;

lens_f = lens_d * F / D;
lens_v = lens_f + v;
micr_N = ceil(D / lens_d);
sen_N_total = micr_N * sen_N;
sen_d = lens_d / sen_N;

if mod(sen_N_total, 2) == 0
    sen_N_total = sen_N_total + 1;
end

im_max = zeros(3,1);

%% ===================== 8. 火焰图送入 LF_sim =====================
tic
for k = 1:3
    fprintf('正在计算第 %d 个通道...\n', k);
    obje = double(obj(:,:,k));

    obj_w = 0.12 * 1000;   % mm
    obj_h = 0.12 * 1000;   % mm

    im = LF_sim(obje, d, obj_w, obj_h, D, N_line, F, v, lens_d, lens_f, sen_N);
    im_max(k) = max(im(:));

    save(sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', d, v, N_line, k), 'im');
end
t = toc;
fprintf('\n模拟光场传输时间 t = %.3f s\n', t);

%% ===================== 9. 组装 5D 光场数据 =====================
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

if ~isequal(size(LF_R), size(LF_G), size(LF_B))
    error('三个通道的光场尺寸不一致');
end

LF = cat(5, LF_R, LF_G, LF_B);

fprintf('\n5D 光场 LF 的维度为: ');
disp(size(LF));
fprintf('LF 的维数 ndims(LF) = %d\n', ndims(LF));

save(sprintf('./dataRGB/LF_5D_d_%d_v_%d_Nline_%d.mat', d, v, N_line), 'LF', '-v7.3');
fprintf('5D 光场 LF 已保存完成！\n');

%% ===================== 10. 显示输入场景 =====================
figure;
imshow(obj, []);
title('输入 LF_sim 的几何图形辐射图');

%% ===================== 11. 显示微透镜阵列传感器图像 =====================
im_RGB = reshape4to2_im(d, v, N_line, micr_N, sen_N, im_max);

figure;
imshow(im_RGB, []);
title('带微透镜阵列的传感器图像');

imwrite(im_RGB, sprintf('./dataRGB/im_d_%d_v_%d_Nline_%d_RGB.jpg', d, v, N_line));

fprintf('\n仿真完成，四层几何图形及微透镜图像已生成。\n');