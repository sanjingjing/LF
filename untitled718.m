clear all;
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

%% =========================================================
%  1. 火焰模型参数
%
%  本程序采用五层径向分层模型：
%
%  1）沿相机光轴 x 方向，将火焰直径 2R 连续划分为5层；
%  2）每层厚度为 2R/5，层与层之间没有空隙；
%  3）每层只使用该层中心位置 x_i 的二维温度分布；
%  4）使用层中心温度代表该层整个厚度范围内的温度；
%  5）每层保留有限厚度对应的辐射发射能力；
%  6）每层使用其中心对应的不同物距调用 LF_sim；
%  7）5个分层光场在线性光场数据层面叠加。
%
%  坐标约定：
%  相机位于火焰 +x 方向，光线沿 +x 方向传播。
%% =========================================================

t1 = 1200;           % 温度增量，K
t2 = 900;            % 基础温度，K
m  = 3;              % 温度场形状参数
n  = 0.9;

R = 0.04;            % 火焰半径，m
Z = 0.40;            % 火焰高度，m

Ny = 500;            % 横向网格数
Nz = 500;            % 轴向网格数

[y, z] = meshgrid( ...
    linspace(-R, R, Ny), ...
    linspace(0, Z, Nz));

%% =========================================================
%  2. 单色辐射参数
%
%  所有参数统一使用 SI 单位。
%% =========================================================

lambda = 0.610e-6;   % 波长，m
c1 = 3.7418e-16;     % 第一辐射常数
c2 = 1.4388e-2;      % 第二辐射常数，m*K

%% =========================================================
%  3. 介质辐射参数
%
%  对应纯吸收火焰：
%  单次散射反照率 omega = 0
%  因此散射系数 ks = 0
%% =========================================================

ke = 50;             % 消光系数，m^-1
ks = 0;              % 散射系数，m^-1
ka = ke - ks;        % 吸收系数，m^-1

if ka < 0
    error('吸收系数 ka 不能为负，请检查 ke 和 ks。');
end

%% =========================================================
%  4. 沿 x 方向连续划分为5层
%
%  火焰光轴方向范围：
%
%       -R <= x <= R
%
%  五层划分：
%
%  S1: [-R,    -3R/5]
%  S2: [-3R/5, -R/5 ]
%  S3: [-R/5,   R/5 ]
%  S4: [ R/5,   3R/5]
%  S5: [ 3R/5,  R   ]
%
%  当 R = 0.04 m 时：
%
%  S1: [-40, -24] mm，中心 -32 mm
%  S2: [-24,  -8] mm，中心 -16 mm
%  S3: [ -8,   8] mm，中心   0 mm
%  S4: [  8,  24] mm，中心  16 mm
%  S5: [ 24,  40] mm，中心  32 mm
%% =========================================================

N_layer = 5;

% 每层的边界
x_edges = linspace(-R, R, N_layer + 1);

% 每层中心位置
x_center = 0.5 .* ...
    (x_edges(1:end-1) + x_edges(2:end));

% 每层厚度
dx_layer = diff(x_edges);

fprintf('火焰沿 x 方向连续划分为 %d 层。\n', N_layer);

fprintf('每层厚度为 %.6f m = %.3f mm。\n\n', ...
    dx_layer(1), dx_layer(1) * 1000);

for i = 1:N_layer

    fprintf(['S%d：左边界 = %+.5f m，' ...
             '右边界 = %+.5f m，中心 = %+.5f m\n'], ...
        i, ...
        x_edges(i), ...
        x_edges(i+1), ...
        x_center(i));
end

%% =========================================================
%  5. 计算每层中心位置的温度和辐射
%
%  第 i 层使用中心坐标 x_i：
%
%       T_i(y,z) = T(x_i,y,z)
%
%  即：
%
%  - 沿层厚方向 x，不再继续区分温度变化；
%  - 每层整个厚度使用中心温度代表；
%  - 但温度仍随 y 和 z 变化；
%  - 所以 T_i 是二维温度分布，不是单一温度值。
%% =========================================================

% 保存五层温度场
T_layer = zeros(Nz, Ny, N_layer, 'single');

% 保存五层出射辐射
I_layer = zeros(Nz, Ny, N_layer, 'single');

% 对每一个 y，圆柱截面朝向 +x 方向的边界为：
%
%       x_front = sqrt(R^2-y^2)
%
% 相机位于 +x 方向。
x_front = sqrt(max(R.^2 - y.^2, 0));

for i = 1:N_layer

    %% 5.1 当前层中心坐标

    xi = x_center(i);

    %% 5.2 当前层中心截面的二维温度分布

    term_i = m .* ( ...
        xi.^2 ./ R.^2 + ...
        y.^2  ./ R.^2 + ...
        z.^2  ./ Z.^2) - n;

    T_i = t1 .* exp(-(term_i.^2)) + t2;

    %% 5.3 根据中心温度计算单色黑体辐射

    Eb_i = c1 .* lambda.^(-5) ./ ...
        expm1(c2 ./ (lambda .* T_i));

    %% 5.4 当前体积层与圆柱火焰的实际交叠厚度
    %
    % 当前层几何范围：
    %
    %       x_edges(i) <= x <= x_edges(i+1)
    %
    % 当前 y 位置圆柱截面的范围：
    %
    %       -x_front <= x <= x_front
    %
    % 两者取交集，得到当前 y 位置的实际火焰层厚度。

    x_intersect_left = max( ...
        x_edges(i), ...
        -x_front);

    x_intersect_right = min( ...
        x_edges(i+1), ...
        x_front);

    dx_effective = max( ...
        x_intersect_right - x_intersect_left, ...
        0);

    %% 5.5 当前有限厚度火焰层的等效发射率
    %
    % 虽然温度只取层中心值，但每层仍具有有限厚度。
    %
    % 纯吸收、均匀介质层的等效发射率：
    %
    %       alpha_i = 1-exp(-ka*dx_effective)

    alpha_i = 1 - exp(-ka .* dx_effective);

    %% 5.6 当前层辐射穿过前方火焰介质的距离
    %
    % 相机位于 +x 方向。
    %
    % 当前层朝向相机的一侧边界为：
    %
    %       x_intersect_right
    %
    % 从该边界到火焰前侧边界的距离为：
    %
    %       L_after = x_front-x_intersect_right

    L_after = max( ...
        x_front - x_intersect_right, ...
        0);

    %% 5.7 当前层到达火焰边界的辐射贡献
    %
    % 当前层自身发射：
    %
    %       Eb_i*alpha_i
    %
    % 经过前方火焰介质后：
    %
    %       I_i = Eb_i*alpha_i*exp(-ke*L_after)

    I_i = Eb_i .* ...
          alpha_i .* ...
          exp(-ke .* L_after);

    % 圆柱火焰外部不应产生辐射
    I_i(dx_effective <= 0) = 0;

    %% 5.8 保存结果

    T_layer(:, :, i) = single(T_i);
    I_layer(:, :, i) = single(I_i);

    fprintf('已完成第 %d 层：中心 x = %+.5f m。\n', ...
        i, xi);
end

%% =========================================================
%  6. 五层辐射叠加
%
%  该结果表示五个体积层到达火焰边界的出射辐射之和：
%
%       I_total = I_1+I_2+I_3+I_4+I_5
%% =========================================================

I_total = sum(double(I_layer), 3);

%% =========================================================
%  7. 显示五层中心温度分布
%% =========================================================

figure('Color', 'white');

for i = 1:N_layer

    subplot(1, N_layer, i);

    imagesc( ...
        linspace(-R, R, Ny), ...
        linspace(0, Z, Nz), ...
        double(T_layer(:, :, i)));

    set(gca, 'YDir', 'normal');

    axis image;
    axis off;

    title(sprintf( ...
        'S_%d\nx_c=%+.0f mm', ...
        i, ...
        x_center(i) * 1000));
end

colormap(jet(256));

sgtitle('五层中心截面的温度分布');

%% =========================================================
%  8. 显示五层辐射图和叠加图
%
%  显示操作不改变后续 LF_sim 的实际输入数据。
%% =========================================================

figure('Color', 'white');

for i = 1:N_layer

    subplot(1, N_layer + 1, i);

    imagesc( ...
        linspace(-R, R, Ny), ...
        linspace(0, Z, Nz), ...
        double(I_layer(:, :, i)));

    set(gca, 'YDir', 'normal');

    axis image;
    axis off;

    title(sprintf( ...
        'S_%d\nx_c=%+.0f mm', ...
        i, ...
        x_center(i) * 1000));
end

subplot(1, N_layer + 1, N_layer + 1);

imagesc( ...
    linspace(-R, R, Ny), ...
    linspace(0, Z, Nz), ...
    I_total);

set(gca, 'YDir', 'normal');

axis image;
axis off;

title('五层叠加');

colormap(gray(256));

sgtitle(sprintf( ...
    '五层中心温度近似下的辐射图，k_e=%.1f m^{-1}', ...
    ke));

%% =========================================================
%  9. 使用统一的全局缩放系数
%
%  注意：
%
%  1）不能对各层分别归一化；
%  2）不能分别做 gamma 变换；
%  3）不能分别转换成 uint8；
%
%  否则会破坏五层之间的真实相对辐射强度。
%
%  这里只使用一个统一缩放系数，把所有层缩放到 LF_sim
%  常用的 0～255 数值范围。
%% =========================================================

global_I_max = max(I_total(:));

if global_I_max <= 0
    error('总辐射强度为零，请检查温度场和介质参数。');
end

LF_input_peak = 255;

global_scale = LF_input_peak / global_I_max;

fprintf('\n总辐射最大值 = %.6e\n', global_I_max);

fprintf('统一 LF_sim 输入缩放系数 = %.6e\n', ...
    global_scale);

%% =========================================================
%  10. 光场相机参数
%% =========================================================

% 火焰中心到主透镜的距离
d_center = 600;       % mm

% 主透镜与传感器参数
v = 16;               % 像距，mm
N_line = 81;          % 主透镜离散化数量
d_m = 18;           

D = 4;                % 主透镜直径，mm
F = 16;               % 主透镜焦距，mm

lens_d = 0.01;        % 微透镜直径，mm
sen_N = 20;           % 每个微透镜后的传感器像素数

lens_f = lens_d * F / D;
lens_v = lens_f + v;  

micr_N = ceil(D / lens_d);

sen_N_total = micr_N * sen_N;
sen_d = lens_d / sen_N; 

if mod(sen_N_total, 2) == 0
    sen_N_total = sen_N_total + 1;
end

% LF_sim 输入物体的物理尺寸
%
% 横向范围：-R～R，因此宽度为 2R
% 轴向范围：0～Z，因此高度为 Z

obj_w = 2 * R * 1000; % 80 mm
obj_h = Z * 1000;     % 400 mm

%% =========================================================
%  11. 计算五层对应的不同物距
%
%  相机位于火焰 +x 方向：
%
%       d_i = d_center-x_i
%
%  x 越大，表示该层越靠近相机，物距越小。
%
%  五层中心：
%
%       [-32,-16,0,16,32] mm
%
%  对应物距：
%
%       [632,616,600,584,568] mm
%% =========================================================

d_layer = d_center - x_center .* 1000;

fprintf('\n各层 LF_sim 物距如下：\n');

for i = 1:N_layer

    fprintf(['S%d：中心位置 = %+7.2f mm，' ...
             '物距 d_i = %.3f mm\n'], ...
        i, ...
        x_center(i) * 1000, ...
        d_layer(i));
end

%% =========================================================
%  12. 每层使用不同物距调用 LF_sim
%
%  LF_layer{i}：
%  第 i 层单独发光时对应的光场。
%
%  LF_total：
%  五层同时发光时的总光场。
%% =========================================================

LF_layer = cell(N_layer, 1);

LF_total = [];

layer_LF_max = zeros(N_layer, 1);

tic;

for i = 1:N_layer

    %% 12.1 当前层二维辐射图
    %
    % 所有层采用同一个全局缩放系数。

    obj_i = double(I_layer(:, :, i)) .* global_scale;

    %% 12.2 当前层物距

    d_i = d_layer(i);

    fprintf('\n正在调用 LF_sim：S%d/%d，d=%.3f mm...\n', ...
        i, ...
        N_layer, ...
        d_i);

    %% 12.3 当前层光场成像

    LF_i = double(LF_sim( ...
        obj_i, ...
        d_i, ...
        obj_w, ...
        obj_h, ...
        D, ...
        N_line, ...
        F, ...
        v, ...
        lens_d, ...
        lens_f, ...
        sen_N));

    %% 12.4 检查不同物距下输出尺寸是否一致

    if isempty(LF_total)

        LF_total = zeros(size(LF_i), 'double');

    elseif ~isequal(size(LF_total), size(LF_i))

        error(['第 %d 层 LF_sim 输出尺寸与前面层不同。' ...
               '请检查 LF_sim 是否会随物距改变输出数组尺寸。'], ...
               i);
    end

    %% 12.5 保存当前层光场

    LF_layer{i} = LF_i;

    %% 12.6 在线性光场数据层面叠加

    LF_total = LF_total + LF_i;

    layer_LF_max(i) = max(LF_i(:));

    %% 12.7 单独保存当前层光场

    save(sprintf( ...
        '.\\dataRGB\\LF_centerTemp_S%d_d_%.3f.mat', ...
        i, ...
        d_i), ...
        'LF_i', ...
        'i', ...
        'd_i', ...
        '-v7.3');
end

elapsed_time = toc;

fprintf('\n五层 LF_sim 计算完成，总耗时 %.3f s。\n', ...
    elapsed_time);

%% =========================================================
%  13. 组装总五维光场数据
%
%  当前为 610 nm 单色辐射模型。
%
%  因此，将同一个单色总光场复制到 R、G、B 三个通道：
%
%       LF(:,:,:,:,1) = LF_total
%       LF(:,:,:,:,2) = LF_total
%       LF(:,:,:,:,3) = LF_total
%
%  若 LF_total 为四维：
%
%       [U,V,S,T]
%
%  则 LF 为五维：
%
%       [U,V,S,T,3]
%% =========================================================

LF = cat(5, ...
    LF_total, ...
    LF_total, ...
    LF_total);

fprintf('\n单色总光场 LF_total 的尺寸为：\n');
disp(size(LF_total));

fprintf('RGB 五维光场 LF 的尺寸为：\n');
disp(size(LF));

fprintf('LF 的维数 ndims(LF) = %d\n', ...
    ndims(LF));

%% =========================================================
%  14. 保存五层模型的全部数据
%% =========================================================

save(sprintf( ...
    '.\\dataRGB\\LF_5layers_centerTemp_dcenter_%d_v_%d_Nline_%d.mat', ...
    d_center, ...
    v, ...
    N_line), ...
    'LF', ...
    'LF_total', ...
    'LF_layer', ...
    'T_layer', ...
    'I_layer', ...
    'I_total', ...
    'x_edges', ...
    'x_center', ...
    'dx_layer', ...
    'd_layer', ...
    'ke', ...
    'ka', ...
    'ks', ...
    'lambda', ...
    'global_scale', ...
    'obj_w', ...
    'obj_h', ...
    '-v7.3');

fprintf('五层中心温度模型的数据已保存。\n');

%% =========================================================
%  15. 保存分层参数表
%% =========================================================

LayerIndex = (1:N_layer).';

XLeft_m = x_edges(1:end-1).';
XRight_m = x_edges(2:end).';
XCenter_m = x_center.';

LayerThickness_m = dx_layer.';
ObjectDistance_mm = d_layer.';

LFMax = layer_LF_max;

layer_info = table( ...
    LayerIndex, ...
    XLeft_m, ...
    XRight_m, ...
    XCenter_m, ...
    LayerThickness_m, ...
    ObjectDistance_mm, ...
    LFMax);

writetable( ...
    layer_info, ...
    '.\dataRGB\layer_info_5layers_centerTemp.csv');

fprintf('分层位置和物距参数表已保存。\n');

%% =========================================================
%  16. 保存兼容 reshape4to2_im 的总光场文件
%
%  原 reshape4to2_im 函数可能会根据以下文件名读取数据：
%
%  im_d_600_v_16_Nline_81_1.mat
%  im_d_600_v_16_Nline_81_2.mat
%  im_d_600_v_16_Nline_81_3.mat
%
%  这里保存的是五层叠加后的总光场。
%
%  文件名中的 d_center=600 只用于兼容原函数，
%  不代表五层全部位于600 mm平面。
%% =========================================================

im = LF_total;

for k = 1:3

    save(sprintf( ...
        '.\\dataRGB\\im_d_%d_v_%d_Nline_%d_%d.mat', ...
        d_center, ...
        v, ...
        N_line, ...
        k), ...
        'im', ...
        '-v7.3');
end

im_max = repmat( ...
    max(LF_total(:)), ...
    3, ...
    1);

%% =========================================================
%  17. 显示五层叠加后的传感器图像
%% =========================================================

im_RGB = reshape4to2_im( ...
    d_center, ...
    v, ...
    N_line, ...
    micr_N, ...
    sen_N, ...
    im_max);

figure('Color', 'white');

imshow(im_RGB, []);

title('五层中心温度模型的总光场传感器图像');

imwrite( ...
    im_RGB, ...
    '.\dataRGB\sensor_5layers_centerTemp_RGB.png');

fprintf('\n已生成五层叠加后的总传感器图像。\n');

disp('The end!');