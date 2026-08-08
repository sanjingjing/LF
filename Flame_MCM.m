clear;
clc;
close all;

addpath(genpath('.\image'));
addpath(genpath('.\image2'));
addpath(genpath('.\data'));
addpath(genpath('.\dataRGB'));
addpath(genpath('.\len_sub'));
addpath(genpath('.\频域重聚焦精简'));

if ~exist('.\dataRGB', 'dir')
    mkdir('.\dataRGB');
end

if ~exist('.\data7layers', 'dir')
    mkdir('.\data7layers');
end

%% ============================================================
%  1. 按论文建立七层火焰模型
%
%  坐标：
%  Xw：光场相机光轴方向，也是火焰分层方向
%  Yw：火焰轴向，即竖直方向
%  Zw：火焰横向
%
%  论文参数：
%  R = 0.04 m
%  Z = 0.40 m
%  七层中心位置：
%  -0.03、-0.02、-0.01、0、0.01、0.02、0.03 m
%  每层厚度：
%  0.01 m
%% ============================================================

T_base = 1200;             % 基础温度，K
T_delta = 600;             % 温度增量，K

shape_m = 3;
shape_n = 0.9;

RF = 0.04;                 % 火焰径向尺寸 R，m
HF = 0.40;                 % 火焰轴向尺寸 Z，m

N_layer = 7;

% 七层中心位置，单位 m
x_layer = -0.03 : 0.01 : 0.03;

% 每层厚度，单位 m
layer_thickness = 0.01;
half_thickness = layer_thickness / 2;

if length(x_layer) ~= N_layer
    error('七层中心位置数量不正确。');
end

% 二维火焰层的离散尺寸
%
% 正式计算：
Ny = 401;
Nz = 161;

% 初次调试时可以暂时使用：
% Ny = 101;
% Nz = 41;

y_vec = linspace(0, HF, Ny);
z_vec = linspace(-RF, RF, Nz);

dy = y_vec(2) - y_vec(1);
dz = z_vec(2) - z_vec(1);

% 每一层是 Y-Z 平面上的二维火焰分布
[Yw, Zw] = ndgrid(y_vec, z_vec);

%% ============================================================
%  2. 七层温度场
%
%  论文公式：
%
%  T(r,z) = 1200 + 600*exp( ...
%      -[3*((r/R)^2+(z/Z)^2)-0.9]^2 )
%
%  在本程序坐标中：
%
%  r^2 = Xw^2 + Zw^2
%  z   = Yw
%
%  每一层厚度范围内，使用该层中心 x_layer(k)
%  对应的温度分布。
%% ============================================================

T_layers = zeros(Ny, Nz, N_layer);

% 当前层与圆柱火焰实际相交的厚度
effective_thickness = zeros(Ny, Nz, N_layer);

% 每层到达相机方向的有效辐射
source_layers = zeros(Ny, Nz, N_layer);

%% ============================================================
%  3. 单色黑体辐射参数
%% ============================================================

lambda_R = 0.610;          % 单色波长，um

c1 = 3.7418e8;             % W·um^4/m^2
c2 = 1.4388e4;             % um·K

% 无散射条件
ka = 30;                   % 吸收系数，m^-1
ke = ka;               
%% ============================================================
%  蒙特卡洛辐射传输参数
%
%  当前仍采用：
%  ke = ka
%
%  即只考虑粒子发射和吸收，不考虑散射。
%% ============================================================

% 固定随机数种子，保证每次运行结果可以重复
rng(20260729, 'twister');

% 单位投影面积上的抽样光线数，单位：ray/m^2
%
% 当前网格：
% dy = 0.001 m
% dz = 0.0005 m
% 单个网格面积约为 5e-7 m^2
%
% ray_density = 4e7 时，
% 每个有效网格约抽样 20 条光线。
ray_density = 4e7;

% 每个Y-Z网格的投影面积
cell_area = dy * dz;

% 每个有效网格中的蒙特卡洛光线数量
N_ray_cell = max( ...
    1, ...
    round(ray_density * cell_area));

fprintf('\n蒙特卡洛参数：\n');
fprintf('单元面积：%.6e m^2\n', cell_area);
fprintf('抽样光线密度：%.6e ray/m^2\n', ray_density);
fprintf('每个有效网格抽样光线数：%d\n', N_ray_cell);
%% ============================================================
%  4. 分别计算七层温度和辐射贡献
%% ============================================================
%  4. 使用蒙特卡洛法计算七层温度和辐射贡献
%
%  每层独立计算：
%
%  1. 使用层中心位置计算该层温度；
%  2. 在层厚范围内均匀随机产生发射位置；
%  3. 随机抽样粒子自由程；
%  4. 判断光线被吸收还是离开火焰；
%  5. 将成功离开火焰的光线能量累加。
%
%  LF_sim 不做任何修改。
%% ============================================================

% 对每一个横向位置 Zw，
% 计算圆柱火焰沿 Xw 方向的半弦长
half_chord = sqrt(max(RF^2 - Zw.^2, 0));

% 相机位于 Xw 负方向，
% 所以火焰朝向相机的前边界为负半弦
x_front = -half_chord;

% 蒙特卡洛结果
source_layers = zeros(Ny, Nz, N_layer);

% 用解析公式计算的参考结果
% 仅用于检查蒙特卡洛计算是否正确
source_layers_reference = zeros(Ny, Nz, N_layer);

% 记录每层的光线逃逸率
MC_escape_ratio = zeros(N_layer, 1);

% 记录每层蒙特卡洛结果和解析结果的相对误差
MC_relative_error = zeros(N_layer, 1);

for k = 1:N_layer

    %% --------------------------------------------------------
    % 1. 当前层的几何范围
    %% --------------------------------------------------------

    xk = x_layer(k);

    % 当前层的名义左右边界
    layer_left = ...
        xk - half_thickness;

    layer_right = ...
        xk + half_thickness;

    % 当前层与圆柱火焰的实际交叠范围
    x_left_effective = max( ...
        layer_left, ...
        -half_chord);

    x_right_effective = min( ...
        layer_right, ...
        half_chord);

    % 当前层在每个Y-Z网格中的实际厚度
    delta_x_effective = max( ...
        x_right_effective - x_left_effective, ...
        0);

    %% --------------------------------------------------------
    % 2. 当前层温度
    %
    % 在该有限厚度层内，
    % 温度统一使用层中心位置 xk 对应的值。
    %% --------------------------------------------------------

    radial_position = ...
        (xk^2 + Zw.^2) ./ RF^2;

    axial_position = ...
        Yw.^2 ./ HF^2;

    T_k = ...
        T_base + ...
        T_delta .* exp( ...
        -( ...
        shape_m .* ...
        (radial_position + axial_position) ...
        - shape_n ...
        ).^2);

    % 当前层有效区域
    layer_mask = ...
        (delta_x_effective > 0) & ...
        (Yw >= 0) & ...
        (Yw <= HF);

    % 火焰区域外不发射
    T_k(~layer_mask) = T_base;

    %% --------------------------------------------------------
    % 3. 当前层黑体单色辐射强度
    %% --------------------------------------------------------

    Eb_k = ...
        c1 .* lambda_R.^(-5) ./ ...
        ( ...
        exp(c2 ./ (lambda_R .* T_k)) ...
        - 1 ...
        );

    Ib_k = Eb_k ./ pi;

    Ib_k(~layer_mask) = 0;

    %% --------------------------------------------------------
    % 4. 蒙特卡洛抽样光线能量
    %
    % 发射系数为：
    %
    % j = ka * Ib
    %
    % 一个网格在整个层厚中的发射能量为：
    %
    % E_cell = ka*Ib*delta_x*dy*dz
    %% --------------------------------------------------------

    emitted_energy_cell = ...
        ka .* ...
        Ib_k .* ...
        delta_x_effective .* ...
        cell_area;

    % 每条蒙特卡洛抽样光线代表的能量
    ray_energy = ...
        emitted_energy_cell ./ N_ray_cell;

    % 当前层成功离开火焰的总能量
    layer_escape_energy = ...
        zeros(Ny, Nz);

    total_valid_rays = ...
        N_ray_cell * nnz(layer_mask);

    total_escape_rays = 0;

    %% --------------------------------------------------------
    % 5. 蒙特卡洛循环
    %% --------------------------------------------------------

    for i_ray = 1:N_ray_cell

        % 在当前有限厚度层内均匀随机产生发射位置
        random_emission_position = ...
            rand(Ny, Nz);

        x_emission = ...
            x_left_effective + ...
            random_emission_position .* ...
            delta_x_effective;

        % 从发射位置到火焰前表面的传播距离
        propagation_distance = ...
            x_emission - x_front;

        propagation_distance = max( ...
            propagation_distance, ...
            0);

        %% 随机抽样粒子自由程

        random_free_path = rand(Ny, Nz);

        % 防止随机数恰好为0
        random_free_path = max( ...
            random_free_path, ...
            eps);

        if ke > 0

            free_path = ...
                -log(random_free_path) ./ ke;

        else

            % 无衰减时，自由程视为无穷大
            free_path = ...
                inf(Ny, Nz);
        end

        %% 判断光线是否成功离开火焰

        escape_mask = ...
            (free_path >= propagation_distance) & ...
            layer_mask;

        % 将成功离开火焰的光线能量累加
        layer_escape_energy = ...
            layer_escape_energy + ...
            ray_energy .* double(escape_mask);

        total_escape_rays = ...
            total_escape_rays + ...
            nnz(escape_mask);
    end

    %% --------------------------------------------------------
    % 6. 从网格能量转换回二维辐射强度
    %
    % LF_sim 接收的是二维辐射图，
    % 因此除以Y-Z网格面积。
    %% --------------------------------------------------------

    layer_emission_MC = ...
        layer_escape_energy ./ cell_area;

    layer_emission_MC(~layer_mask) = 0;

    %% --------------------------------------------------------
    % 7. 解析结果，仅用于校验蒙特卡洛模型
    %% --------------------------------------------------------

    % 从火焰前表面到当前层前边界的距离
    foreground_length = max( ...
        x_left_effective - x_front, ...
        0);

    transmission_reference = ...
        exp(-ka .* foreground_length);

    layer_emission_reference = ...
        Ib_k .* ...
        (1 - exp( ...
        -ka .* delta_x_effective)) .* ...
        transmission_reference;

    layer_emission_reference(~layer_mask) = 0;

    %% --------------------------------------------------------
    % 8. 保存当前层结果
    %% --------------------------------------------------------

    T_layers(:, :, k) = T_k;

    effective_thickness(:, :, k) = ...
        delta_x_effective;

    % 后续进入LF_sim的是蒙特卡洛结果
    source_layers(:, :, k) = ...
        layer_emission_MC;

    source_layers_reference(:, :, k) = ...
        layer_emission_reference;

    %% --------------------------------------------------------
    % 9. 蒙特卡洛误差和逃逸率
    %% --------------------------------------------------------

    if total_valid_rays > 0

        MC_escape_ratio(k) = ...
            total_escape_rays / total_valid_rays;
    end

    reference_norm = ...
        norm(layer_emission_reference(:));

    if reference_norm > 0

        MC_relative_error(k) = ...
            norm( ...
            layer_emission_MC(:) - ...
            layer_emission_reference(:)) ./ ...
            reference_norm;
    end

    fprintf('\n第%d层蒙特卡洛计算完成：\n', k);

    fprintf( ...
        '层中心：%.3f m\n', ...
        xk);

    fprintf( ...
        '名义范围：[%.3f, %.3f] m\n', ...
        layer_left, ...
        layer_right);

    fprintf( ...
        '物距：%.3f mm\n', ...
        800 + xk * 1000);

    fprintf( ...
        '抽样光线数：%d\n', ...
        total_valid_rays);

    fprintf( ...
        '成功逃逸光线数：%d\n', ...
        total_escape_rays);

    fprintf( ...
        '光线逃逸率：%.4f\n', ...
        MC_escape_ratio(k));

    fprintf( ...
        '与解析结果的相对误差：%.4f\n', ...
        MC_relative_error(k));
end
%% ============================================================
%% ============================================================
%  蒙特卡洛结果与解析结果对比
%% ============================================================

% 七层蒙特卡洛辐射叠加
MC_projection = ...
    sum(source_layers, 3);

% 七层解析辐射叠加
reference_projection = ...
    sum(source_layers_reference, 3);

%% 归一化

MC_projection_norm = ...
    MC_projection - min(MC_projection(:));

if max(MC_projection_norm(:)) > 0

    MC_projection_norm = ...
        MC_projection_norm ./ ...
        max(MC_projection_norm(:));
end

reference_projection_norm = ...
    reference_projection - ...
    min(reference_projection(:));

if max(reference_projection_norm(:)) > 0

    reference_projection_norm = ...
        reference_projection_norm ./ ...
        max(reference_projection_norm(:));
end

%% MSE

MC_difference = ...
    MC_projection_norm - ...
    reference_projection_norm;

MSE_MC = ...
    mean(MC_difference(:).^2);

fprintf('\n蒙特卡洛七层总图评价：\n');
fprintf('MSE = %.8f\n', MSE_MC);

%% SSIM，需要Image Processing Toolbox

if exist('ssim', 'file') == 2

    SSIM_MC = ssim( ...
        MC_projection_norm, ...
        reference_projection_norm);

    fprintf('SSIM = %.8f\n', SSIM_MC);

else

    SSIM_MC = NaN;

    fprintf( ...
        '当前环境没有ssim函数，跳过SSIM计算。\n');
end

%% 显示对比结果

figure('Color', 'white', ...
    'Name', '蒙特卡洛结果与解析结果对比');

subplot(1, 3, 1);

imagesc( ...
    z_vec, ...
    y_vec, ...
    reference_projection_norm);

set(gca, 'YDir', 'normal');

axis image;
axis tight;

title('解析辐射结果');

xlabel('Z_w (m)');
ylabel('Y_w (m)');

colorbar;

subplot(1, 3, 2);

imagesc( ...
    z_vec, ...
    y_vec, ...
    MC_projection_norm);

set(gca, 'YDir', 'normal');

axis image;
axis tight;

title('蒙特卡洛辐射结果');

xlabel('Z_w (m)');
ylabel('Y_w (m)');

colorbar;

subplot(1, 3, 3);

imagesc( ...
    z_vec, ...
    y_vec, ...
    abs(MC_difference));

set(gca, 'YDir', 'normal');

axis image;
axis tight;

title('绝对误差');

xlabel('Z_w (m)');
ylabel('Y_w (m)');

colorbar;

colormap(jet(256));

%% 保存蒙特卡洛计算结果

save('.\data7layers\seven_layer_MC_source.mat', ...
    'source_layers', ...
    'source_layers_reference', ...
    'T_layers', ...
    'effective_thickness', ...
    'MC_escape_ratio', ...
    'MC_relative_error', ...
    'MSE_MC', ...
    'SSIM_MC', ...
    'ray_density', ...
    'N_ray_cell', ...
    'x_layer', ...
    'layer_thickness', ...
    '-v7.3');
%  5. 显示七层温度分布
%% ============================================================

figure('Color', 'white', ...
    'Name', '七层火焰温度分布');

for k = 1:N_layer

    subplot(2, 4, k);

    T_show = T_layers(:, :, k);

    layer_mask_show = ...
        effective_thickness(:, :, k) > 0;

    T_show(~layer_mask_show) = NaN;

    imagesc(z_vec, y_vec, T_show);

    set(gca, 'YDir', 'normal');

    axis image;
    axis tight;

    xlabel('横向 Z_w (m)');
    ylabel('轴向 Y_w (m)');

    title(sprintf( ...
        '第%d层，x=%.3f m', ...
        k, x_layer(k)));

    colorbar;
end

colormap(jet(256));

%% ============================================================
%  6. 显示七层有效辐射分布
%% ============================================================

figure('Color', 'white', ...
    'Name', '七层火焰有效辐射分布');

source_max = max(source_layers(:));

for k = 1:N_layer

    subplot(2, 4, k);

    layer_show = source_layers(:, :, k);

    if source_max > 0
        layer_show = layer_show ./ source_max;
    end

    imagesc(z_vec, y_vec, layer_show);

    set(gca, 'YDir', 'normal');

    axis image;
    axis tight;

    xlabel('横向 Z_w (m)');
    ylabel('轴向 Y_w (m)');

    title(sprintf( ...
        '第%d层，x=%.3f m', ...
        k, x_layer(k)));

    colorbar;
end

colormap(gray(256));

%% ============================================================
%  7. 七层叠加后的二维辐射图
%
%  该图只用来检查七层火焰辐射模型。
%  光场仿真时仍然逐层调用 LF_sim。
%% ============================================================

I_projection_raw = sum(source_layers, 3);

I_projection = ...
    I_projection_raw - min(I_projection_raw(:));

if max(I_projection(:)) > 0
    I_projection = ...
        I_projection ./ max(I_projection(:));
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

title('七层火焰辐射叠加图');

cb = colorbar;
cb.Label.String = '归一化单色辐射强度';

%% ============================================================
%  8. 光场相机参数
%% ============================================================

% 论文中火焰中心线距离主透镜 0.8 m
d_center = 800;            % mm

v = 16;                    % 主透镜像距，mm

% 正式仿真
N_line = 81;

% 初次调试可以暂时改成
% N_line = 21;

D = 4;                     % 主透镜直径，mm
F = 16;                    % 主透镜焦距，mm

lens_d = 0.01;             % 微透镜直径，mm
sen_N = 20;                % 每个微透镜后像素数

lens_f = lens_d * F / D;

% 与 LF_sim 保持一致
micr_N = fix(D / lens_d);

% 每个二维火焰层的物理尺寸
obj_w = 2 * RF * 1000;     % 80 mm
obj_h = HF * 1000;         % 400 mm

fprintf('\n二维火焰层物理尺寸：\n');
fprintf('obj_w = %.2f mm\n', obj_w);
fprintf('obj_h = %.2f mm\n', obj_h);

fprintf('\n七层物距：\n');

for k = 1:N_layer
    fprintf('第%d层：%.3f mm\n', ...
        k, d_center + x_layer(k) * 1000);
end

%% ============================================================
%  9. 七层分别通过光场相机成像
%
%  重要：
%  每层禁止单独归一化，否则会破坏层间能量比例。
%% ============================================================

im_volume = [];

valid_layer_number = 0;

layer_LF_max = zeros(N_layer, 1);

% 是否保存每一层单独的四维光场
save_single_layer_LF = true;

tic;

for k = 1:N_layer

    layer_image = ...
        double(source_layers(:, :, k));

    if max(layer_image(:)) <= 0

        warning('第%d层没有有效辐射，跳过。', k);

        continue;
    end

    % 当前层实际物距，单位 mm
    d_slice = ...
        d_center + x_layer(k) * 1000;

    fprintf('\n正在计算七层火焰第 %d/%d 层\n', ...
        k, N_layer);

    fprintf('层中心位置：%.3f m\n', x_layer(k));
    fprintf('当前物距：%.3f mm\n', d_slice);

    %% 当前层四维光场

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

    %% 初始化完整七层光场

    if isempty(im_volume)

        im_volume = ...
            zeros(size(im_layer), 'single');
    end

    if ~isequal(size(im_volume), size(im_layer))

        error( ...
            '第%d层 LF_sim 输出尺寸与前面层不一致。', ...
            k);
    end

    %% 记录当前层最大值

    layer_LF_max(k) = max(im_layer(:));

    %% 保存当前单层四维光场

   if save_single_layer_LF

    im = single(im_layer); 

    save_name = fullfile( ...
        'data7layers', ...
        sprintf('LF_layer_%02d.mat', k));

    save( ...
        save_name, ...
        'im', ...
        'd_slice', ...
        'layer_image', ...
        '-v7.3');
   end
    %% 七层光场按能量线性叠加

    im_volume = ...
        im_volume + single(im_layer);

    valid_layer_number = ...
        valid_layer_number + 1;

    clear im_layer;
end

simulation_time = toc;

fprintf('\n有效火焰层数：%d\n', ...
    valid_layer_number);

fprintf('七层火焰光场计算时间：%.3f s\n', ...
    simulation_time);

if isempty(im_volume)

    error( ...
        '没有获得有效光场，请检查 source_layers 或 LF_sim。');
end

fprintf('\n七层叠加四维光场尺寸：\n');
disp(size(im_volume));

%% 保存完整七层四维光场

save_name = fullfile( ...
    'data7layers', ...
    'LF_7layer_sum.mat');

save( ...
    save_name, ...
    'im_volume', ...
    'x_layer', ...
    'layer_thickness', ...
    'd_center', ...
    'layer_LF_max', ...
    '-v7.3');

%% ============================================================
%  10. 按原程序格式保存 RGB 三个通道
%
%  当前是 0.610 um 单色模型，
%  因此三个通道暂时使用同一个四维光场。
%% ============================================================

im_max = zeros(3, 1);

for channel = 1:3

    im = im_volume; 
    im_max(channel) = max(im(:));

   outputDir = fullfile('.', 'dataRGB');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

filename = sprintf( ...
    'im_d_%d_v_%d_Nline_%d_%d.mat', ...
    d_center, ...
    v, ...
    N_line, ...
    channel);

save(fullfile(outputDir, filename), 'im', '-v7.3');

    fprintf( ...
        '第%d通道保存完成，最大值 %.6e\n', ...
        channel, ...
        im_max(channel));
end

%% ============================================================
%  11. 可选：组装五维 RGB 光场
%  五维数据占用内存很大。
%  内存充足时将 save_5D_LF 设置为 true。

file_name = sprintf( ...
    'LF_4D_7layer_d_%d_v_%d_Nline_%d.mat', ...
    d_center, ...
    v, ...
    N_line);

file_path = fullfile('dataRGB', file_name);

save( ...
    file_path, ...
    'im_volume', ...
    '-v7.3');
%% ============================================================
%  12. 显示七层火焰辐射叠加图
obj_gray = ...
    uint8(I_projection_display * 255);

obj_RGB = ...
    repmat(obj_gray, [1, 1, 3]);

figure('Color', 'white');

imshow(obj_RGB);

title('七层火焰沿光轴方向的辐射叠加图');

%% ============================================================
%  13. 生成微透镜阵列传感器图像

im_RGB = reshape4to2_im( ...
    d_center, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure('Color', 'white');

imshow(im_RGB);

title('七层火焰的微透镜阵列传感器图像');

outputDir = fullfile('.', 'dataRGB');

% 文件夹不存在时自动创建
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fileName = sprintf( ...
    'im_7layer_d_%d_v_%d_Nline_%d_RGB.png', ...
    d_center, ...
    v, ...
    N_line);

outputPath = fullfile(outputDir, fileName);

imwrite(im_RGB, outputPath);

fprintf('图像已保存至：%s\n', outputPath);
fprintf('\n七层火焰光场成像仿真完成。\n');

disp('The end!');