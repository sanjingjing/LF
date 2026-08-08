clear ;
clc;
close all;

hsize = 37;
num_iter = 50;
image_size= 512;
imageFolder = fullfile(pwd, 'flame_layers_512');

%% ===================== 1. 火焰辐射模型 =====================
Q1 = imread(fullfile( imageFolder, 'layer_4_x_+0.00_m.png'));

Q2 = imread(fullfile(imageFolder, 'layer_2_x_-0.02_m.png'));
Q3 = imread(fullfile(imageFolder,  'layer_6_x_+0.02_m.png'));

%% 转换为double类型
Q1 = im2double(Q1);
Q2 = im2double(Q2);
Q3 = im2double(Q3);

%% 三张图片叠加
I = Q1 + Q2 + Q3;
I = mat2gray(I);
figure;
imshow(I, []);
title('叠加后的原始图像');


%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.29);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.81);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 9.73);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.36);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.17);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.59);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 9.05);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 6.19);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator(size(I), hsize, 2.61);

%% Forward imaging model

I1 = Op11(Q1) + Op12(Q2) + Op13(Q3);
I2 = Op21(Q1) + Op22(Q2) + Op23(Q3);
I3 = Op31(Q1) + Op32(Q2) + Op33(Q3);

%% Stack measured images b = [b1; b2; b3]

image = zeros(1536, 512);
image(1:512, :)       = I1;
image(513:1024, :)    = I2;
image(1025:1536, :)   = I3;

image = max(real(image), 0);

figure;
imshow(mat2gray(image));
title('观测叠加图像 b');

%% Define block forward operator A

A = @(x) [
    Op11(x(1:512, :)) + Op12(x(513:1024, :)) + Op13(x(1025:1536, :));
    Op21(x(1:512, :)) + Op22(x(513:1024, :)) + Op23(x(1025:1536, :));
    Op31(x(1:512, :)) + Op32(x(513:1024, :)) + Op33(x(1025:1536, :))
];

%% Define adjoint operator A_adj

A_adj = @(y) [
    Op11_adj(y(1:512, :)) + Op21_adj(y(513:1024, :)) + Op31_adj(y(1025:1536, :));
    Op12_adj(y(1:512, :)) + Op22_adj(y(513:1024, :)) + Op32_adj(y(1025:1536, :));
    Op13_adj(y(1:512, :)) + Op23_adj(y(513:1024, :)) + Op33_adj(y(1025:1536, :))
];

%% Initial estimate

init = ones(size(image)) .* mean(image(:));

%% Richardson-Lucy deconvolution

field_rl = deconvRL_operator(A, A_adj, image, num_iter, init);

field_rl = real(field_rl);
field_rl(field_rl < 0) = 0;


figure;
imshow(field_rl, []);
title('RL去卷积整体结果');

%% Split reconstructed layers
rec1 = field_rl(1:512, :);
rec2 = field_rl(513:1024, :);
rec3 = field_rl(1025:1536, :);

figure;
subplot(3,3,1); imshow(mat2gray(Q1)); title('原始图像1');
subplot(3,3,2); imshow(mat2gray(Q2)); title('原始图像2');
subplot(3,3,3); imshow(mat2gray(Q3)); title('原始图像3');

subplot(3,3,4);imshow(I1, []);title('重聚焦层1');
subplot(3,3,5);imshow(I2, []);title('重聚焦层2');
subplot(3,3,6);imshow(I3, []);title('重聚焦层3');
subplot(3,3,7);imshow(rec1, []);title('恢复层 x1');
subplot(3,3,8);imshow(rec2, []);title('恢复层 x2');
subplot(3,3,9);imshow(rec3, []);title('恢复层 x3');
%% ===================== MATLAB官方SSIM计算 =====================

% 确保输入数据为double类型
Q1 = double(Q1);
Q2 = double(Q2);
Q3 = double(Q3);

rec1 = double(rec1);
rec2 = double(rec2);
rec3 = double(rec3);

% 原始图像由im2double得到，理论动态范围为[0,1]
dynamicRange = 1;

% 火焰区域阈值：
% 原始图像强度大于该层最大值的1%时，判定为火焰区域
thresholdRatio = 0.01;

%% --------------------- 第1层SSIM ---------------------

% globalSSIM1为全图平均SSIM
% ssimMap1为局部SSIM分布图
[globalSSIM1, ssimMap1] = ssim(rec1, Q1, ...
    'DynamicRange', dynamicRange, ...
    'Radius', 1.5, ...
    'Exponents', [1 1 1]);

% 根据原始图像Q1提取火焰区域
maxQ1 = max(Q1(:));

if maxQ1 <= 0
    error('Q1中没有有效火焰像素。');
end

flameMask1 = Q1 > thresholdRatio * maxQ1;

% 防止部分MATLAB版本返回的SSIM图尺寸与原图不同
if ~isequal(size(flameMask1), size(ssimMap1))
    flameMask1 = imresize(flameMask1, ...
        [size(ssimMap1, 1), size(ssimMap1, 2)], 'nearest');
end

% 取出火焰区域内有效的局部SSIM值
validSSIM1 = ssimMap1(flameMask1);
validSSIM1 = validSSIM1(isfinite(validSSIM1));

if isempty(validSSIM1)
    error('Q1对应的火焰区域内没有有效SSIM值。');
end

% 只计算火焰区域平均SSIM
flameSSIM1 = mean(validSSIM1);


%% --------------------- 第2层SSIM ---------------------

[globalSSIM2, ssimMap2] = ssim(rec2, Q2, ...
    'DynamicRange', dynamicRange, ...
    'Radius', 1.5, ...
    'Exponents', [1 1 1]);

maxQ2 = max(Q2(:));

if maxQ2 <= 0
    error('Q2中没有有效火焰像素。');
end

flameMask2 = Q2 > thresholdRatio * maxQ2;

if ~isequal(size(flameMask2), size(ssimMap2))
    flameMask2 = imresize(flameMask2, ...
        [size(ssimMap2, 1), size(ssimMap2, 2)], 'nearest');
end

validSSIM2 = ssimMap2(flameMask2);
validSSIM2 = validSSIM2(isfinite(validSSIM2));

if isempty(validSSIM2)
    error('Q2对应的火焰区域内没有有效SSIM值。');
end

flameSSIM2 = mean(validSSIM2);


%% --------------------- 第3层SSIM ---------------------

[globalSSIM3, ssimMap3] = ssim(rec3, Q3, ...
    'DynamicRange', dynamicRange, ...
    'Radius', 1.5, ...
    'Exponents', [1 1 1]);

maxQ3 = max(Q3(:));

if maxQ3 <= 0
    error('Q3中没有有效火焰像素。');
end

flameMask3 = Q3 > thresholdRatio * maxQ3;

if ~isequal(size(flameMask3), size(ssimMap3))
    flameMask3 = imresize(flameMask3, ...
        [size(ssimMap3, 1), size(ssimMap3, 2)], 'nearest');
end

validSSIM3 = ssimMap3(flameMask3);
validSSIM3 = validSSIM3(isfinite(validSSIM3));

if isempty(validSSIM3)
    error('Q3对应的火焰区域内没有有效SSIM值。');
end

flameSSIM3 = mean(validSSIM3);


%% ===================== 输出SSIM计算结果 =====================

fprintf('\n================ SSIM计算结果 ================\n');

fprintf('Q1 与 rec1：\n');
fprintf('  全图平均SSIM   = %.6f\n', globalSSIM1);
fprintf('  火焰区域SSIM   = %.6f\n\n', flameSSIM1);

fprintf('Q2 与 rec2：\n');
fprintf('  全图平均SSIM   = %.6f\n', globalSSIM2);
fprintf('  火焰区域SSIM   = %.6f\n\n', flameSSIM2);

fprintf('Q3 与 rec3：\n');
fprintf('  全图平均SSIM   = %.6f\n', globalSSIM3);
fprintf('  火焰区域SSIM   = %.6f\n', flameSSIM3);

fprintf('==============================================\n');


%% ===================== 使用imagesc绘制SSIM Map =====================

figure('Color', 'w', ...
       'Name', '三层局部SSIM分布图', ...
       'Position', [100 100 1300 430]);

% 第1层SSIM Map
subplot(1, 3, 1);

imagesc(ssimMap1, [-1 1]);
axis image;
axis off;
colorbar;

title(sprintf(['Q1 与 rec1 的局部SSIM图\n', ...
               '火焰区域平均SSIM = %.4f'], ...
               flameSSIM1), ...
      'FontSize', 11);

% 第2层SSIM Map
subplot(1, 3, 2);

imagesc(ssimMap2, [-1 1]);
axis image;
axis off;
colorbar;

title(sprintf(['Q2 与 rec2 的局部SSIM图\n', ...
               '火焰区域平均SSIM = %.4f'], ...
               flameSSIM2), ...
      'FontSize', 11);

% 第3层SSIM Map
subplot(1, 3, 3);

imagesc(ssimMap3, [-1 1]);
axis image;
axis off;
colorbar;

title(sprintf(['Q3 与 rec3 的局部SSIM图\n', ...
               '火焰区域平均SSIM = %.4f'], ...
               flameSSIM3), ...
      'FontSize', 11);

sgtitle('原始火焰层与对应RL恢复层的局部SSIM分布', ...
        'FontSize', 14, ...
        'FontWeight', 'bold');


%% ===================== 火焰区域平均SSIM柱状图 =====================

% 柱状图只使用火焰区域的平均SSIM
flameSSIMValues = [flameSSIM1, flameSSIM2, flameSSIM3];

figure('Color', 'w', ...
       'Name', '火焰区域平均SSIM柱状图', ...
       'Position', [200 150 700 520]);

bar(1:3, flameSSIMValues, 0.60);

set(gca, ...
    'XTick', 1:3, ...
    'XTickLabel', {'Q1-rec1', 'Q2-rec2', 'Q3-rec3'}, ...
    'FontSize', 12, ...
    'LineWidth', 1);

xlabel('原始火焰层与对应恢复层', ...
       'FontSize', 13);

ylabel('火焰区域平均SSIM', ...
       'FontSize', 13);

title('三层RL恢复结果的火焰区域SSIM', ...
      'FontSize', 14);

grid on;
box on;

% 根据SSIM是否存在负值设置坐标范围
if min(flameSSIMValues) >= 0
    ylim([0 1.08]);
else
    ylim([-1 1.08]);
end

% 在柱顶标注具体数值
text(1, flameSSIM1 + 0.025, ...
    sprintf('%.4f', flameSSIM1), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

text(2, flameSSIM2 + 0.025, ...
    sprintf('%.4f', flameSSIM2), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');

text(3, flameSSIM3 + 0.025, ...
    sprintf('%.4f', flameSSIM3), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 11, ...
    'FontWeight', 'bold');