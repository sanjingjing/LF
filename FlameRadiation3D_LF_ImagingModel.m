clear;
clc;
close all;

addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\频域重聚焦精简'));

%% ============================================================
%  1. 三维圆柱形火焰物理模型
%  坐标定义：
%  Xw：相机光轴方向，即火焰深度方向
%  Yw：火焰轴向，即竖直方向
%  Zw：火焰横向
%% ============================================================

T_delta = 1200;       % 温度增量，K
T_ambient = 900;      % 环境温度，K
shape_m = 3;
shape_n = 0.9;

RF = 0.02;            % 火焰半径，m
HF = 0.20;            % 火焰轴向高度，m

% 三维体素数
% 调试时建议 Nx=21；正式计算可增加至41、61
Nx = 21;              % 深度方向，每层厚度维2RF/(Nx-1)
Ny = 401;             % 火焰高度方向
Nz = 161;             % 横向方向

x_vec = linspace(-RF, RF, Nx);
y_vec = linspace(0, HF, Ny);
z_vec = linspace(-RF, RF, Nz);

dx = x_vec(2) - x_vec(1);
dy = y_vec(2) - y_vec(1);
dz = z_vec(2) - z_vec(1);

% ndgrid输出顺序分别对应Y、Z、X
[Yw, Zw, Xw] = ndgrid(y_vec, z_vec, x_vec);

% 圆柱形火焰区域
flame_mask = ...
    (Xw.^2 + Zw.^2 <= RF^2) & ...
    (Yw >= 0) & (Yw <= HF);

% 论文式（5-1）三维温度分布
normalized_position = ...
    (Xw.^2 + Zw.^2) ./ RF^2 + ...
    Yw.^2 ./ HF^2;

T = T_delta .* exp( ...
    -(shape_m .* normalized_position - shape_n).^2 ...
    ) + T_ambient;

% 火焰外部不作为参与性介质
T(~flame_mask) = T_ambient;

fprintf('三维温度场尺寸：\n');
disp(size(T));

fprintf('三维温度范围：%.2f K ～ %.2f K\n', ...
    min(T(:)), max(T(:)));
%  2. 显示三维火焰温度场
%% ============================================================

% 当前T的维度为：
% T(y, z, x) = [Ny, Nz, Nx]
%
% MATLAB的slice/interp3需要：
% V(y, x, z) = [Ny, Nx, Nz]
%
% 因此交换第2维和第3维
T_plot = permute(T, [1, 3, 2]);

fprintf('用于slice显示的温度场尺寸：\n');
disp(size(T_plot));
% 应输出：[Ny, Nx, Nz]，即401 × 41 × 161

figure('Color', 'white');

% 使用坐标向量形式，避免NDGRID和MESHGRID结构冲突
slice( ...
    x_vec, ...                         % X：火焰深度方向
    y_vec, ...                         % Y：火焰轴向方向
    z_vec, ...                         % Z：火焰横向方向
    T_plot, ...
    [-0.75*RF, 0, 0.75*RF], ...       % X方向切面
    [0.04, 0.10, 0.16], ...           % Y方向切面
    0);                               % Z方向切面

shading interp;
axis equal;
axis tight;
view(3);
grid on;

xlabel('深度坐标 X_w (m)');
ylabel('轴向坐标 Y_w (m)');
zlabel('横向坐标 Z_w (m)');
title('三维火焰温度场');

cb = colorbar;
cb.Label.String = '温度 T (K)';

colormap(jet(256));

%% ============================================================
%  3. 单色黑体辐射模型
%% ============================================================

lambda_R = 0.610;     % 单色波长，um

% 当波长单位采用um时，c1必须使用下面的单位
c1 = 3.7418e8;        % W·um^4/m^2
c2 = 1.4388e4;        % um·K

% 黑体单色辐射出射度
Eb = c1 .* lambda_R.^(-5) ./ ...
    (exp(c2 ./ (lambda_R .* T)) - 1);

% 转换为黑体单色辐射强度
Ib = Eb ./ pi;

% 火焰外部没有发射
Ib(~flame_mask) = 0;

fprintf('单色黑体辐射强度范围：%.6e ～ %.6e\n', ...
    min(Ib(:)), max(Ib(:)));

%% ============================================================
%  4. 吸收与衰减模型
%
% 这里采用无散射辐射传输方程：
%
% dI/ds = ka*Ib - ke*I
%
% 当前假定：
% ke = ka，即忽略散射，仅考虑吸收和自身发射。
%% ============================================================

ke = 30;               % 衰减系数，m^-1
ka = ke;               % 无散射条件下，吸收系数等于衰减系数

% 每一个深度体素内的吸收系数
ka_volume = ka .* double(flame_mask);

% 每层光线传播到火焰前表面时的透射率
transmission_to_camera = ones(Ny, Nz);

% 存放各深度层到达火焰前表面的辐射强度
source_layers = zeros(Ny, Nz, Nx);

% 假定相机位于Xw负方向
% 因此从x=-RF开始，依次走向x=+RF
for ix = 1:Nx

    ka_slice = ka_volume(:, :, ix);
    Ib_slice = Ib(:, :, ix);

    % 当前体素层自身向相机方向贡献的辐射强度
    % 采用体素内辐射传输解析解：
    % dI = Ib*(1-exp(-ka*dx))*前方透射率
    layer_emission = ...
        Ib_slice .* ...
        (1 - exp(-ka_slice .* dx)) .* ...
        transmission_to_camera;

    source_layers(:, :, ix) = layer_emission;

    % 更新之后各层到相机的透射率
    transmission_to_camera = ...
        transmission_to_camera .* exp(-ka_slice .* dx);
end

%% ============================================================
%  5. 显示三维火焰沿光轴积分后的二维辐射图
%
% 注意：
% 该图只用于检查三维辐射模型是否正确。
% 后面不会直接把它作为唯一平面送入LF_sim，
% 而是把每个深度层分别送入LF_sim。
%% ============================================================

I_projection = sum(source_layers, 3);

I_projection = I_projection - min(I_projection(:));

if max(I_projection(:)) > 0
    I_projection = I_projection ./ max(I_projection(:));
end

I_projection_display = I_projection .^ 0.7;

figure('Color', 'white');

imagesc(z_vec, y_vec, I_projection_display);
set(gca, 'YDir', 'normal');

axis image;
axis tight;
colormap(gray(256));

xlabel('横向坐标 Z_w (m)');
ylabel('轴向坐标 Y_w (m)');
title('三维火焰辐射传输积分图');

cb = colorbar;
cb.Label.String = '归一化单色辐射强度';

%% ============================================================
%  6. 光场相机参数
%% ============================================================

% 主透镜到火焰中心的距离
d_center = 600;        % mm

v = 16;                % 主透镜像距，mm
N_line = 81;           % 主透镜离散化数

D = 4;                 % 主透镜直径，mm
F = 16;                % 主透镜焦距，mm

lens_d = 0.01;         % 微透镜直径，mm
sen_N = 20;            % 每个微透镜对应的像素数

lens_f = lens_d * F / D;
lens_v = lens_f + v;

micr_N = ceil(D / lens_d);
sen_N_total = micr_N * sen_N;
sen_d = lens_d / sen_N;

if mod(sen_N_total, 2) == 0
    sen_N_total = sen_N_total + 1;
end

% 每个二维深度层的物理尺寸
obj_w = 2 * RF * 1000;    % 横向宽度，mm
obj_h = HF * 1000;        % 纵向高度，mm

fprintf('\n二维深度层物理尺寸：\n');
fprintf('横向宽度 obj_w = %.2f mm\n', obj_w);
fprintf('轴向高度 obj_h = %.2f mm\n', obj_h);

%% ============================================================
%  7. 将三维火焰各个深度层分别送入LF_sim
%
% 关键修改：
% 原程序只有一个固定x截面；
% 这里每个x位置都具有不同物距d_slice。
%% ============================================================

im_volume = [];
valid_layer_number = 0;

tic;

for ix = 1:Nx

    layer_image = source_layers(:, :, ix);

    % 忽略没有明显辐射贡献的截面
    if max(layer_image(:)) <= 0
        continue;
    end

    % 每一层单独归一化会破坏层间能量比例，因此只能统一缩放
    layer_image = double(layer_image);

    % 当前火焰层到主透镜的实际距离
    % Xw<0为靠近相机的一侧
    d_slice = d_center + x_vec(ix) * 1000;

    fprintf('正在计算三维火焰第 %d/%d 层，物距 %.3f mm\n', ...
        ix, Nx, d_slice);

    % 将当前真实深度层送入光场相机
    im_layer = LF_sim( ...
        layer_image, ...
        d_slice, ...
        obj_w, ...
        obj_h, ...
        D, ...
        N_line, ...
        F, ...
        v, ...
        lens_d, ...
        lens_f, ...
        sen_N);

    % 将不同深度产生的光场能量进行线性叠加
    if isempty(im_volume)
        im_volume = zeros(size(im_layer), 'double');
    end

    if ~isequal(size(im_volume), size(im_layer))
        error('不同深度层LF_sim输出尺寸不一致。');
    end

    im_volume = im_volume + double(im_layer);
    valid_layer_number = valid_layer_number + 1;
end

simulation_time = toc;

fprintf('\n三维火焰有效深度层数：%d\n', valid_layer_number);
fprintf('三维火焰光场传播时间：%.3f s\n', simulation_time);

if isempty(im_volume)
    error('没有得到有效光场结果，请检查source_layers或LF_sim。');
end

%% ============================================================
%  8. 保存三维火焰四维光场
%
% reshape4to2_im要求：
% 文件名：
% im_d_d_v_v_Nline_Nline_k.mat
%
% 变量名：
% im
%% ============================================================

if isempty(im_volume)
    error('im_volume为空，三维火焰光场计算失败。');
end

% 检查LF_sim输出是否为四维
if ndims(im_volume) ~= 4
    error('LF_sim累加结果必须是四维数组，当前尺寸为：%s', ...
        mat2str(size(im_volume)));
end

fprintf('\n三维火焰累加光场尺寸：\n');
disp(size(im_volume));

% 不要在这里转换为uint8
% 保存double类型，reshape4to2_im中统一归一化
im_max = zeros(3, 1);

for k = 1:3

    % 当前为0.610 um单色模型，RGB三个通道暂时相同
    im = im_volume;

    im_max(k) = max(im(:));

    save(sprintf( ...
        './dataRGB/im_d_%d_v_%d_Nline_%d_%d.mat', ...
        d_center, v, N_line, k), ...
        'im', '-v7.3');

    fprintf('第%d通道光场数据保存完成，最大值为 %.6e\n', ...
        k, im_max(k));
end

%% ============================================================
%  9. 组装5D光场
%% ============================================================

LF_R = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_1.mat', ...
    d_center, v, N_line), 'im');

LF_G = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_2.mat', ...
    d_center, v, N_line), 'im');

LF_B = load(sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_3.mat', ...
    d_center, v, N_line), 'im');

LF_R = LF_R.im;
LF_G = LF_G.im;
LF_B = LF_B.im;

if ~isequal(size(LF_R), size(LF_G), size(LF_B))
    error('RGB三个通道的四维光场尺寸不一致。');
end

LF = cat(5, LF_R, LF_G, LF_B);

fprintf('\n5D光场LF尺寸为：\n');
disp(size(LF));

save(sprintf( ...
    './dataRGB/LF_5D_d_%d_v_%d_Nline_%d.mat', ...
    d_center, v, N_line), ...
    'LF', '-v7.3');

fprintf('三维火焰5D光场保存完成。\n');

%% ============================================================
%  10. 显示三维火焰辐射积分图
%% ============================================================

obj_gray = uint8(I_projection_display * 255);
obj_RGB = repmat(obj_gray, [1, 1, 3]);

figure('Color', 'white');
imshow(obj_RGB);
title('三维火焰沿光轴方向的辐射积分图');

%% ============================================================
%  11. 使用原来的reshape4to2_im生成传感器图像
%% ============================================================

im_RGB = reshape4to2_im( ...
    d_center, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure('Color', 'white');
imshow(im_RGB);
title('三维火焰的微透镜阵列传感器图像');

imwrite(im_RGB, sprintf( ...
    './dataRGB/im_d_%d_v_%d_Nline_%d_RGB.png', ...
    d_center, v, N_line));

fprintf('\n三维火焰传感器图像绘制完成。\n');
disp('The end!');