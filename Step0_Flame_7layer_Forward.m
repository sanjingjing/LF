clear;
clc;
close all;

%% ============================================================
%  七层火焰直接叠加光场仿真
%
%  本版本已经删除：
%  1. 蒙特卡洛随机抽样；
%  2. 自由程和逃逸率计算；
%  3. 蒙特卡洛与解析结果对比；
%  4. MSE、SSIM 等蒙特卡洛评价。
%
%  现在采用：
%  每层直接辐射贡献 = ka * Ib * 有效层厚
%  然后七层分别通过 LF_sim 成像，并进行线性叠加。
%% ============================================================

%% 0. 添加代码路径并创建输出文件夹

codeFolders = { ...
    'image', ...
    'image2', ...
    'len_sub', ...
    '频域重聚焦精简'};

for iFolder = 1:numel(codeFolders)
    folderPath = fullfile('.', codeFolders{iFolder});

    if exist(folderPath, 'dir')
        addpath(genpath(folderPath));
    end
end

outputDirRGB = fullfile('.', 'dataRGB');
outputDir7 = fullfile('.', 'data7layers');

if ~exist(outputDirRGB, 'dir')
    mkdir(outputDirRGB);
end

if ~exist(outputDir7, 'dir')
    mkdir(outputDir7);
end

%% ============================================================
%  1. 建立七层火焰模型
%
%  坐标：
%  Xw：光场相机光轴方向，也是火焰分层方向
%  Yw：火焰轴向，即竖直方向
%  Zw：火焰横向
%% ============================================================

T_base = 1200;             % 基础温度，K
T_delta = 600;             % 温度增量，K

shape_m = 3;
shape_n = 0.9;

RF = 0.04;                 % 火焰径向尺寸，m
HF = 0.40;                 % 火焰轴向尺寸，m

N_layer = 7;

% 七层中心位置，单位 m
x_layer = -0.03 : 0.01 : 0.03;

% 每层厚度，单位 m
layer_thickness = 0.01;
half_thickness = layer_thickness / 2;

if numel(x_layer) ~= N_layer
    error('七层中心位置数量不正确。');
end

% 二维火焰层离散尺寸
Ny = 401;
Nz = 161;

% 初次调试可改小：
% Ny = 101;
% Nz = 41;

y_vec = linspace(0, HF, Ny);
z_vec = linspace(-RF, RF, Nz);

[Yw, Zw] = ndgrid(y_vec, z_vec);

%% ============================================================
%  2. 单色黑体辐射参数
%% ============================================================

lambda_R = 0.610;          % 单色波长，um

c1 = 3.7418e8;             % W·um^4/m^2
c2 = 1.4388e4;             % um·K

ka = 30;                   % 吸收/发射系数，m^-1

%% ============================================================
%  3. 直接计算七层温度场和辐射贡献
%
%  不再使用蒙特卡洛方法。
%
%  当前采用光学薄层近似：
%
%  source = ka * Ib * delta_x
%
%  其中：
%  Ib      为单色黑体辐射强度；
%  delta_x 为当前层与圆柱火焰相交的实际厚度。
%
%  七层之间不计算前景吸收，直接线性叠加。
%% ============================================================

T_layers = zeros(Ny, Nz, N_layer);
effective_thickness = zeros(Ny, Nz, N_layer);
source_layers = zeros(Ny, Nz, N_layer);

% 圆柱火焰在每个 Zw 位置沿 Xw 方向的半弦长
half_chord = sqrt(max(RF^2 - Zw.^2, 0));

for k = 1:N_layer

    xk = x_layer(k);

    %% 当前层的名义范围

    layer_left = xk - half_thickness;
    layer_right = xk + half_thickness;

    %% 当前层与圆柱火焰的实际交叠范围

    x_left_effective = max(layer_left, -half_chord);
    x_right_effective = min(layer_right, half_chord);

    delta_x_effective = max( ...
        x_right_effective - x_left_effective, ...
        0);

    %% 当前层温度场

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

    layer_mask = ...
        (delta_x_effective > 0) & ...
        (Yw >= 0) & ...
        (Yw <= HF);

    % 火焰区域外不参与辐射
    T_k(~layer_mask) = T_base;

    %% 当前层单色黑体辐射强度

    Eb_k = ...
        c1 .* lambda_R.^(-5) ./ ...
        ( ...
        exp(c2 ./ (lambda_R .* T_k)) ...
        - 1 ...
        );

    Ib_k = Eb_k ./ pi;
    Ib_k(~layer_mask) = 0;

    %% 当前层直接辐射贡献

    layer_source = ...
        ka .* ...
        Ib_k .* ...
        delta_x_effective;

    layer_source(~layer_mask) = 0;

    %% 保存当前层结果

    T_layers(:, :, k) = T_k;
    effective_thickness(:, :, k) = delta_x_effective;
    source_layers(:, :, k) = layer_source;

    fprintf('\n第 %d 层直接辐射计算完成：\n', k);
    fprintf('层中心：%.3f m\n', xk);
    fprintf('名义范围：[%.3f, %.3f] m\n', ...
        layer_left, layer_right);
    fprintf('最大辐射值：%.6e\n', ...
        max(layer_source(:)));
end

%% 保存七层直接辐射结果

sourceFile = fullfile( ...
    outputDir7, ...
    'seven_layer_direct_source.mat');

save( ...
    sourceFile, ...
    'source_layers', ...
    'T_layers', ...
    'effective_thickness', ...
    'x_layer', ...
    'layer_thickness', ...
    'ka', ...
    'lambda_R', ...
    '-v7.3');

fprintf('\n七层直接辐射数据已保存：%s\n', sourceFile);

%% ============================================================
%  4. 显示七层温度分布
%% ============================================================

figure( ...
    'Color', 'white', ...
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
%  5. 显示七层直接辐射分布
%% ============================================================

figure( ...
    'Color', 'white', ...
    'Name', '七层火焰直接辐射分布');

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
%  6. 七层物面辐射直接叠加
%% ============================================================

I_projection_raw = sum(source_layers, 3);

I_projection = ...
    I_projection_raw - min(I_projection_raw(:));

if max(I_projection(:)) > 0
    I_projection = ...
        I_projection ./ max(I_projection(:));
end

% 仅用于显示的伽马调整
I_projection_display = I_projection .^ 0.7;

figure( ...
    'Color', 'white', ...
    'Name', '七层火焰辐射直接叠加图');

imagesc(z_vec, y_vec, I_projection_display);

set(gca, 'YDir', 'normal');

axis image;
axis tight;

colormap(gray(256));

xlabel('横向坐标 Z_w (m)');
ylabel('轴向坐标 Y_w (m)');

title('七层火焰辐射直接叠加图');

cb = colorbar;
cb.Label.String = '归一化单色辐射强度';

%% ============================================================
%  7. 光场相机参数
%% ============================================================

d_center = 800;            % 火焰中心物距，mm
v = 16;                    % 主透镜像距，mm

N_line = 81;               % 正式仿真
% N_line = 21;             % 初次调试

D = 4;                     % 主透镜直径，mm
F = 16;                    % 主透镜焦距，mm

lens_d = 0.02;             % 微透镜直径，mm
sen_N = 20;                % 每个微透镜后像素数

lens_f = lens_d * F / D;
micr_N = fix(D / lens_d);

obj_w = 2 * RF * 1000;     % 80 mm
obj_h = HF * 1000;         % 400 mm

fprintf('\n二维火焰层物理尺寸：\n');
fprintf('obj_w = %.2f mm\n', obj_w);
fprintf('obj_h = %.2f mm\n', obj_h);

fprintf('\n七层物距：\n');

for k = 1:N_layer
    fprintf( ...
        '第%d层：%.3f mm\n', ...
        k, ...
        d_center + x_layer(k) * 1000);
end

%% ============================================================

%  8. 七层分别通过 LF_sim 成像，并在四维光场层面叠加
%
%  每一层先使用自己的物距调用 LF_sim，得到一个四维光场；
%  然后执行：
%
%      im_sum = im_sum + im_layer;
%
%  得到七层叠加后的四维光场。
%
%  注意：每层不能单独归一化，否则会破坏层间能量比例。
%% ============================================================

im_sum = [];
valid_layer_number = 0;
layer_LF_max = zeros(N_layer, 1);

% 是否保存每一层单独的四维光场
save_single_layer_LF = true;

tic;

for k = 1:N_layer

    layer_image = double(source_layers(:, :, k));

    if max(layer_image(:)) <= 0
        warning('第%d层没有有效辐射，跳过。', k);
        continue;
    end

    % 当前层自己的物距，单位 mm
    d_slice = ...
        d_center + x_layer(k) * 1000;

    fprintf('\n正在计算第 %d/%d 层光场\n', ...
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

    % 转成 single，降低叠加和保存时的内存占用
    im_layer = single(im_layer);

    %% 初始化七层叠加四维光场

    if isempty(im_sum)
        im_sum = zeros(size(im_layer), 'single');
    end

    if ~isequal(size(im_sum), size(im_layer))
        error( ...
            '第%d层 LF_sim 输出尺寸与前面层不一致。', ...
            k);
    end

    layer_LF_max(k) = max(im_layer(:));

    %% 可选：保存当前单层四维光场

    if save_single_layer_LF

        im = im_layer;

        singleLayerFile = fullfile( ...
            outputDir7, ...
            sprintf('LF_layer_%02d.mat', k));

        save( ...
            singleLayerFile, ...
            'im', ...
            'd_slice', ...
            'layer_image', ...
            '-v7.3');
    end

    %% 七层四维光场线性叠加

    im_sum = im_sum + im_layer;

    valid_layer_number = ...
        valid_layer_number + 1;

    clear im_layer;
end

simulation_time = toc;

fprintf('\n有效火焰层数：%d\n', ...
    valid_layer_number);

fprintf('七层火焰光场计算时间：%.3f s\n', ...
    simulation_time);

if isempty(im_sum)
    error( ...
        '没有获得有效光场，请检查 source_layers 或 LF_sim。');
end

fprintf('\n七层叠加后四维光场 im_sum 的尺寸：\n');
disp(size(im_sum));

fprintf('im_sum 的维数 ndims(im_sum) = %d\n', ...
    ndims(im_sum));

%% 保存七层叠加后的四维光场

sumLFFile = fullfile( ...
    outputDir7, ...
    'LF_7layer_4D_sum.mat');

save( ...
    sumLFFile, ...
    'im_sum', ...
    'x_layer', ...
    'layer_thickness', ...
    'd_center', ...
    'layer_LF_max', ...
    '-v7.3');

fprintf('七层叠加四维光场已保存：%s\n', ...
    sumLFFile);

%% ============================================================
%  9. 组装并保存五维 RGB 光场数组 LF
%
%  五维数组的前四维是光场的空间/角度维度，
%  第五维是颜色通道：
%
%      LF(:,:,:,:,1) = R 通道四维光场
%      LF(:,:,:,:,2) = G 通道四维光场
%      LF(:,:,:,:,3) = B 通道四维光场
%
%  当前程序使用 0.610 um 单色火焰模型，因此暂时令
%  R、G、B 三个通道都等于同一个七层叠加四维光场。
%% ============================================================

im_total_R = im_sum;
im_total_G = im_sum;
im_total_B = im_sum;

if ~isequal( ...
        size(im_total_R), ...
        size(im_total_G), ...
        size(im_total_B))

    error( ...
        '三个通道的四维光场尺寸不一致，无法组装五维 LF。');
end

% 沿第 5 维拼接三个四维光场
LF = cat( ...
    5, ...
    im_total_R, ...
    im_total_G, ...
    im_total_B);

fprintf('\n五维光场数组 LF 的尺寸：\n');
disp(size(LF));

fprintf('LF 的维数 ndims(LF) = %d\n', ...
    ndims(LF));

if ndims(LF) ~= 5 || size(LF, 5) ~= 3
    error( ...
        'LF 组装失败：期望得到第五维长度为3的五维数组。');
end

%% 保存五维光场数组，变量名必须是 LF

LF5DFile = fullfile( ...
    outputDirRGB, ...
    sprintf( ...
    'LF_5D_7layer_d_%d_v_%d_Nline_%d.mat', ...
    d_center, ...
    v, ...
    N_line));

save( ...
    LF5DFile, ...
    'LF', ...
    'x_layer', ...
    'layer_thickness', ...
    'd_center', ...
    'v', ...
    'N_line', ...
    '-v7.3');

fprintf('五维光场数组 LF 已保存：%s\n', ...
    LF5DFile);

%% ============================================================
%  10. 按 reshape4to2_im 所需格式保存三个四维通道
%
%  reshape4to2_im 会按照固定文件名读取变量 im，
%  所以这里仍然分别保存 R、G、B 三个四维光场。
%% ============================================================

im_max = zeros(3, 1);

for channel = 1:3

    switch channel
        case 1
            im = im_total_R;
        case 2
            im = im_total_G;
        case 3
            im = im_total_B;
    end

    im_max(channel) = max(im(:));

    filename = sprintf( ...
        'im_d_%d_v_%d_Nline_%d_%d.mat', ...
        d_center, ...
        v, ...
        N_line, ...
        channel);

    channelFile = fullfile( ...
        outputDirRGB, ...
        filename);

    save( ...
        channelFile, ...
        'im', ...
        '-v7.3');

    fprintf( ...
        '第%d通道四维光场保存完成，最大值 %.6e\n', ...
        channel, ...
        im_max(channel));
end

%% ============================================================
%  11. 显示七层火焰物面叠加图


obj_gray = uint8(I_projection_display * 255);
obj_RGB = repmat(obj_gray, [1, 1, 3]);

figure( ...
    'Color', 'white', ...
    'Name', '七层火焰物面叠加图');

imshow(obj_RGB);

title('七层火焰沿光轴方向的直接辐射叠加图');

%% ============================================================
%  12. 生成微透镜阵列传感器图像
%% ============================================================

im_RGB = reshape4to2_im( ...
    d_center, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure( ...
    'Color', 'white', ...
    'Name', '微透镜阵列传感器图像');

imshow(im_RGB);

title('七层火焰的微透镜阵列传感器图像');

fileName = sprintf( ...
    'im_7layer_direct_d_%d_v_%d_Nline_%d_RGB.png', ...
    d_center, ...
    v, ...
    N_line);

outputPath = fullfile( ...
    outputDirRGB, ...
    fileName);

imwrite(im_RGB, outputPath);

fprintf('图像已保存至：%s\n', outputPath);
fprintf('\n七层火焰四维叠加及五维 LF 保存完成。\n');

disp('The end!');
