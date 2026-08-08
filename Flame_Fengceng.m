clear;
close all;
clc;

%% ===================== 0. 显示及保存开关 =====================

showAxes  = true;      % true：显示横纵坐标；false：隐藏坐标轴
saveImages = true;     % true：保存7张512×512灰度图
targetSize = 512;      % 最终图像尺寸

%% ===================== 1. 参数设置 =====================

t1 = 600;              % 温度增量，K
t2 = 1200;             % 环境温度，K
m  = 3;                % 形状参数
n  = 0.9;

R = 0.04;              % 火焰径向尺度，m
Z = 0.40;              % 火焰轴向尺度，m

% 原始 y-z 平面离散
Ny = 500;
Nz = 500;

yVector = linspace(-0.04, 0.04, Ny);
zVector = linspace(0, 0.40, Nz);

[y, z] = meshgrid(yVector, zVector);

%% ===================== 2. 七层切片参数 =====================

% 七层中心坐标
xCenters = -0.03 : 0.01 : 0.03;

numberOfLayers = length(xCenters);

% 每层中心左右各0.005 m，总厚度0.01 m
halfThickness = 0.005;

% 每层厚度方向采样数
NxLayer = 21;

% 保存原始计算结果
I_layers = zeros(Nz, Ny, numberOfLayers);
T_layers = zeros(Nz, Ny, numberOfLayers);

%% ===================== 3. 单色辐射参数 =====================

lambda = 0.610;        % 波长，um

c1 = 3.7418e8;         % W·um^4/m^2
c2 = 1.4388e4;         % um·K

%% ===================== 4. 计算七层火焰 =====================

for layerIndex = 1:numberOfLayers

    xCenter = xCenters(layerIndex);

    xMin = xCenter - halfThickness;
    xMax = xCenter + halfThickness;

    xLocal = linspace(xMin, xMax, NxLayer);

    T_stack  = zeros(Nz, Ny, NxLayer);
    Eb_stack = zeros(Nz, Ny, NxLayer);

    for xIndex = 1:NxLayer

        xCurrent = xLocal(xIndex);

        % 火焰温度分布
        term = m .* ( ...
            xCurrent.^2 ./ R.^2 + ...
            y.^2        ./ R.^2 + ...
            z.^2        ./ Z.^2) - n;

        T = t1 .* exp(-(term.^2)) + t2;

        % 单色黑体辐射模型
        Eb = c1 .* lambda.^(-5) ./ ...
            (exp(c2 ./ (lambda .* T)) - 1);

        T_stack(:, :, xIndex)  = T;
        Eb_stack(:, :, xIndex) = Eb;

    end

    % 当前层厚度范围内的平均温度
    T_layers(:, :, layerIndex) = ...
        trapz(xLocal, T_stack, 3) ./ (xMax - xMin);

    % 当前层厚度范围内的单色辐射积分
    % 不使用Beer衰减
    I_layers(:, :, layerIndex) = ...
        trapz(xLocal, Eb_stack, 3);

end

%% ===================== 5. 构造512×512物理等比例画布 =====================

% z方向物理长度为0.4 m。
% 为形成正方形物理画布，y方向也设置为0.4 m。
zPhysicalLength = zVector(end) - zVector(1);

yCanvasMin = -zPhysicalLength / 2;
yCanvasMax =  zPhysicalLength / 2;

yVector512 = linspace(yCanvasMin, yCanvasMax, targetSize);
zVector512 = linspace(zVector(1), zVector(end), targetSize);

[yCanvas, zCanvas] = meshgrid(yVector512, zVector512);

% 512×512×7结果，初始值为0，即黑色
I_layers_512 = zeros(targetSize, targetSize, numberOfLayers);
T_layers_512 = zeros(targetSize, targetSize, numberOfLayers);

for layerIndex = 1:numberOfLayers

    % 辐射场插值到512×512物理坐标画布
    % 超出原火焰计算范围的区域赋值为0
    I_layers_512(:, :, layerIndex) = interp2( ...
        yVector, ...
        zVector, ...
        I_layers(:, :, layerIndex), ...
        yCanvas, ...
        zCanvas, ...
        'linear', ...
        0);

    % 温度场插值
    T_layers_512(:, :, layerIndex) = interp2( ...
        yVector, ...
        zVector, ...
        T_layers(:, :, layerIndex), ...
        yCanvas, ...
        zCanvas, ...
        'linear', ...
        0);

end

fprintf('原始辐射矩阵大小：%d × %d\n', ...
    size(I_layers, 1), size(I_layers, 2));

fprintf('扩充后辐射矩阵大小：%d × %d\n', ...
    size(I_layers_512, 1), size(I_layers_512, 2));

%% ===================== 6. 统一颜色范围 =====================

% 0必须对应黑色填充区域
radiationMin = 0;
radiationMax = max(I_layers_512(:));

if radiationMax <= 0
    error('辐射强度最大值小于或等于0，请检查参数。');
end

%% ===================== 7. 居中显示七层辐射切片 =====================

figureHandle = figure( ...
    'Color', 'white', ...
    'Units', 'pixels', ...
    'Position', [100, 100, 1100, 670]);

movegui(figureHandle, 'center');

% 2×8布局，每幅图占两列
layout = tiledlayout(2, 8, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

axesHandles = gobjects(numberOfLayers, 1);

% 第一行四张
topTilePositions = [1, 3, 5, 7];

% 第二行三张，居中排列
bottomTilePositions = [10, 12, 14];

tilePositions = [topTilePositions, bottomTilePositions];

for layerIndex = 1:numberOfLayers

    axesHandles(layerIndex) = nexttile( ...
        layout, ...
        tilePositions(layerIndex), ...
        [1, 2]);

    imagesc( ...
        yVector512, ...
        zVector512, ...
        I_layers_512(:, :, layerIndex));

    set(gca, ...
        'YDir', 'normal', ...
        'FontSize', 9);

    % y、z物理范围均为0.4 m，因此显示区域为正方形
    axis image;
    axis tight;
    grid off;

    clim([radiationMin, radiationMax]);

    xCenter = xCenters(layerIndex);
    xMin = xCenter - halfThickness;
    xMax = xCenter + halfThickness;

    title({ ...
        sprintf('第 %d 层', layerIndex), ...
        sprintf('x_c = %.2f m', xCenter), ...
        sprintf('x \\in [%.3f, %.3f] m', xMin, xMax) ...
        }, ...
        'FontSize', 10, ...
        'FontWeight', 'normal');

    %% 横纵坐标显示开关
    if showAxes

        axis on;

        xlabel('径向坐标 y (m)', ...
            'FontSize', 10);

        ylabel('轴向坐标 z (m)', ...
            'FontSize', 10);

        % 可以根据需要修改坐标刻度
        xticks(-0.2 : 0.1 : 0.2);
        yticks(0 : 0.1 : 0.4);

    else

        % 隐藏坐标轴、刻度及坐标名称
        axis off;

    end

end

colormap(gray(256));

% 公共颜色条
cb = colorbar(axesHandles(end));
cb.Layout.Tile = 'east';
cb.Label.String = '厚度积分单色辐射强度';
cb.Label.FontSize = 11;
cb.Label.Rotation = 90;

title(layout, ...
    '七层火焰单色辐射切片（无 Beer 衰减）', ...
    'Interpreter', 'none', ...
    'FontSize', 15, ...
    'FontWeight', 'bold');

%% ===================== 8. 保存512×512原始灰度图 =====================

if saveImages

    outputFolder = 'flame_layers_512';

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    for layerIndex = 1:numberOfLayers

        currentImage = I_layers_512(:, :, layerIndex);

        % 七层共用同一个最大值，保证亮度可以直接比较
        normalizedImage = currentImage ./ radiationMax;

        % 限制在[0,1]
        normalizedImage = max(0, min(1, normalizedImage));
       normalizedImage=flipud(normalizedImage);
        outputFile = fullfile( ...
            outputFolder, ...
            sprintf('layer_%d_x_%+.2f_m.png', ...
            layerIndex, ...
            xCenters(layerIndex)));

        % 直接保存数据矩阵，不包含标题、坐标轴和空白边框
        imwrite(normalizedImage, outputFile);

        fprintf('已保存：%s，尺寸为%d × %d\n', ...
            outputFile, ...
            size(normalizedImage, 1), ...
            size(normalizedImage, 2));

    end

end