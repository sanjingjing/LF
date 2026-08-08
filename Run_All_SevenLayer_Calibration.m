%% 七层火焰前向成像 + 重聚焦参数标定：一键运行
%
% 运行时间可能较长：
% 1. 先生成七层单层四维光场；
% 2. 再对七层分别扫描 alpha 并计算 Brenner 清晰度。

clear;
clc;
close all;

root_folder = fileparts(mfilename('fullpath'));

if isempty(root_folder)
    root_folder = pwd;
end

fprintf('\n========== Step 0：七层火焰光场前向仿真 ==========\n');

run(fullfile( ...
    root_folder, ...
    'Step0_Flame_7layer_Forward.m'));

% Step0 中包含 clear，因此重新获取根目录
root_folder = fileparts(mfilename('fullpath'));

if isempty(root_folder)
    root_folder = pwd;
end

fprintf('\n========== Step 1：七层 alpha 参数标定 ==========\n');

run(fullfile( ...
    root_folder, ...
    'Step1_SevenLayer_Refocus_Brenner_Calibration.m'));

fprintf('\n========== 全部流程完成 ==========\n');
