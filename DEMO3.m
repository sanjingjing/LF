clear all;
clc;
close all;

hsize = 37;
num_iter = 150;

n = 512;

Q1 = imresize(im2double(im2gray(imread('2.png'))), [n n]);
Q2 = imresize(im2double(im2gray(imread('cameraman.tif'))), [n n]);





%% Create Gaussian blur operators

[Op11, Op11_adj, ~] = createGaussianBlurringOperator([n n], hsize, 2.36);
[Op12, Op12_adj, ~] = createGaussianBlurringOperator([n n], hsize, 5.55);
[Op13, Op13_adj, ~] = createGaussianBlurringOperator([n n], hsize, 12.97);

[Op21, Op21_adj, ~] = createGaussianBlurringOperator([n n], hsize, 6.37);
[Op22, Op22_adj, ~] = createGaussianBlurringOperator([n n], hsize, 2.29);
[Op23, Op23_adj, ~] = createGaussianBlurringOperator([n n], hsize, 8.44);

[Op31, Op31_adj, ~] = createGaussianBlurringOperator([n n], hsize, 9.72);
[Op32, Op32_adj, ~] = createGaussianBlurringOperator([n n], hsize, 5.25);
[Op33, Op33_adj, ~] = createGaussianBlurringOperator([n n], hsize, 2.27);

%% Forward imaging model
var = 0.0; % noise variance
I1 = Op11(Q1) + Op12(Q2)+ var * randn([n n])  ;
I2 = Op21(Q1) + Op22(Q2) + var * randn([n n]) ;

figure;
subplot(1,2,1); imshow(mat2gray(I1)); title('叠加图像1');
subplot(1,2,2); imshow(mat2gray(I2)); title('叠加图像2');


%% Stack measured images b = [b1; b2; b3]

image = zeros(1024, 512);
image(1:512, :)       = I1;
image(513:1024, :)    = I2;


image = max(real(image), 0);



%% Define block forward operator A

A = @(x) [
    Op11(x(1:512, :)) + Op12(x(513:1024, :));
    Op21(x(1:512, :)) + Op22(x(513:1024, :) );
    
];

%% Define adjoint operator A_adj

A_adj = @(y) [
    Op11_adj(y(1:512, :)) + Op21_adj(y(513:1024, :));
    Op12_adj(y(1:512, :)) + Op22_adj(y(513:1024, :) );
   
];

%% Initial estimate

init = ones(size(image)) .* mean(image(:));

%% Richardson-Lucy deconvolution

field_rl = deconvRL_operator(A, A_adj, image, num_iter, init);

field_rl = real(field_rl);
field_rl(field_rl < 0) = 0;
field_rl = field_rl ./ max(field_rl(:));



%% Split reconstructed layers

rec1 = field_rl(1:512, :);
rec2 = field_rl(513:1024, :);


figure;
subplot(1,2,1); imshow(rec1, []); title('恢复文字层rec1');
subplot(1,2,2); imshow(rec2, []); title('恢复人像层rec2');
figure;
subplot(1,2,1),imshow(Q1),title('文字');
subplot(1,2,2),imshow(Q2),title('cameraman');
%% ===================== SSIM计算 =====================

K = [0.01 0.03];
L = 1;

% SSIM窗口
window = fspecial('gaussian',[11 11],1.5);


%% -------- 第一层 rec1 与 Q1 --------

[mssim1, ssim_map1] = ssim_github( ...
    Q1, rec1, K, window, L);


%% -------- 第二层 rec2 与 Q2 --------

[mssim2, ssim_map2] = ssim_github( ...
    rec2, Q2, K, window, L);



fprintf('================ SSIM结果 ================\n');

fprintf('rec1 与 Q1 平均SSIM = %.6f\n',mssim1);

fprintf('rec2 与 Q2 平均SSIM = %.6f\n',mssim2);

fprintf('==========================================\n');



%% ===================== SSIM Map显示 =====================

figure('Color','w');


subplot(1,2,1)

imagesc(ssim_map1,[0 1]);
axis image;
axis off;
colorbar;

title(sprintf('rec1 与 Q1 局部SSIM\n平均SSIM=%.4f',mssim1));



subplot(1,2,2)

imagesc(ssim_map2,[0 1]);
axis image;
axis off;
colorbar;

title(sprintf('rec2 与 Q2 局部SSIM\n平均SSIM=%.4f',mssim2));


sgtitle('恢复层与原始层的局部SSIM分布');



%% ===================== SSIM柱状图 =====================

SSIM_value=[mssim1,mssim2];


figure('Color','w');

bar(SSIM_value,0.5);


set(gca,...
    'XTick',1:2,...
    'XTickLabel',{'rec1-Q1','rec2-Q2'},...
    'FontSize',12);


ylabel('平均SSIM');

xlabel('恢复层');

title('RL恢复结果SSIM比较');


ylim([0 1]);

grid on;


% 添加数值

text(1,mssim1+0.03,...
    sprintf('%.4f',mssim1),...
    'HorizontalAlignment','center',...
    'FontWeight','bold');


text(2,mssim2+0.03,...
    sprintf('%.4f',mssim2),...
    'HorizontalAlignment','center',...
    'FontWeight','bold');


