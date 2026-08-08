clear; clc; close all;

%% ========================================================================
%  复刻代码 —— 第六版v2（当前最优配置）
%  对应论文第 4.2 节：光场层析重建火焰三维温度场模拟
%
%  已解决问题(1-11):
%    1: 评价指标 MSE → SSIM          2: 主图辐射场 → 温度场
%    3: 可视化参数 m=2.5 正确        4: 滑动窗口重建
%    5: PSF σ 非线性模型             6: SNR定义(式4-21)
%    7: 普朗克定律 ε/(4π)           8: 积分离散化等价
%    9: 重建中无 Δy_rec             10: 温度域SSIM
%   11: SSIM评价区域→火焰ROI
%
%  第六版修正(12a,12b):
%   12a: 成像+重建PSF → 线性插值（论文原文"线性插值"）
%   12b: Landweber → Van Cittert + β=1/N_win
%        线性PSF使H对称→Van Cittert合法, β=1/N_win→DC一步精确收敛
%
%  未解决(12c): SSIM迭代变化幅度 Δ≈0.04 vs 论文 Δ≈0.40
%  ========================================================================

%% 1. 模拟条件与网格设置 (对应 4.2.1)
Nx = 80; Nz = 120;
X_vec = linspace(-25, 25, Nx);   % 径向 x (mm)
Z_vec = linspace(0, 75, Nz);     % 轴向高度 z (mm)
[X, Z_grid] = meshgrid(X_vec, Z_vec);

% 分层设置：模拟成像用 400 层，重建反演用 7 层
N_sim = 400; y_sim = linspace(-25, 25, N_sim); dy_sim = 50 / (N_sim - 1);
N_rec = 7;   y_rec = [-21, -14, -7, 0, 7, 14, 21]; dy_rec = 7;  % 论文图4-5(b)

% 火焰物理参数 (式 4-17)
t1 = 1200; t2 = 900; m_param = 10; n_param = 0.4; Z_max = 75; R = 25;

% 辐射传输参数 (式 4-18, 4-19)
lambda = 700e-9;          % 700nm波段
c1 = 3.7418e-16;          % 第一辐射常数 (W.m^2)
c2 = 1.4388e-2;           % 第二辐射常数 (m.K)
kappa = 8;                % 吸收系数 m^-1
epsilon = 1;              % soot发射率

% --- PSF 关键参数 ---
sigma_in = 1.0;
sigma_neighbor = 2.8;

% --- 【问题12a】PSF sigma 线性插值模型 ---
% 论文4.2.1节: "其他断层的σ值根据式 4-6 和式 4-9 通过线性插值求解"
% 锚定点: σ(0)=σ_in=1.0, σ(dy_rec)=σ_neighbor=2.8
calc_sigma = @(dy_abs) sigma_in + (sigma_neighbor - sigma_in) * dy_abs / dy_rec;

fprintf('PSF模型: 线性插值\n');
fprintf('  σ(0)=%.1f, σ(dy_rec=%.1fmm)=%.1f, σ(2dy_rec)=%.1f, σ(3dy_rec)=%.1f\n', ...
    calc_sigma(0), dy_rec, calc_sigma(dy_rec), calc_sigma(2*dy_rec), calc_sigma(3*dy_rec));

%% 2. 构建 400 层超采样物理辐射真值场
fprintf('正在计算 400 层的三维温度场与绝对光谱辐射场...\n');
f_sim = zeros(Nz, Nx, N_sim);
for j = 1:N_sim
    y = y_sim(j);
    T = t1 * exp( - ( m_param * ( (Z_grid/Z_max).^2 + (X.^2 + y^2)/(R^2) ) - n_param ).^2 ) + t2;
    I_lambda = (c1 * lambda^-5) ./ (exp(c2 ./ (lambda * T)) - 1) * epsilon / (4 * pi);
    L_m = (25 - y) / 1000;
    f_sim(:,:,j) = I_lambda .* exp(-kappa * L_m);
end

I_max = max(f_sim(:));
f_sim = f_sim / I_max;

%% 3. 前向投影：生成 7 幅光场重聚焦图像 (对应式 4-20)
% ======================================================================
% 成像模型（400层精细积分）：E = Σ_{i=1}^{400} f_i * h_i · Δy_sim
% 这里乘 Δy_sim 因为 f_sim 是辐射强度密度，需要积分得到总辐射。
% 重建模型（3层纯求和）：E = Σ_{j=1}^{N_win} f_j * h_j （无Δy）
% 两者的 f 物理含义不同：成像用密度，重建用"层总贡献"。
% 【问题12a】各子层PSF使用线性插值σ
% ======================================================================
fprintf('正在利用 400 层数据进行前向积分模拟成像...\n');

E_captured = zeros(Nz, Nx, N_rec);
for k = 1:N_rec
    y_focus = y_rec(k);
    E_k = zeros(Nz, Nx);
    for j = 1:N_sim
        dy_abs = abs(y_focus - y_sim(j));
        sigma = calc_sigma(dy_abs);     % 线性插值PSF
        psf_size = ceil(sigma * 3) * 2 + 1;
        h = fspecial('gaussian', [psf_size, psf_size], sigma);
        E_k = E_k + imfilter(f_sim(:,:,j), h, 'replicate', 'same') * dy_sim;
    end
    E_captured(:,:,k) = E_k;
end

% 按论文式 4-21 定义的 SNR 添加高斯白噪声
SNR_dB = 30;
signal_abs_sum = sum(abs(E_captured(:)));
noise_abs_sum  = signal_abs_sum / (10^(SNR_dB/10));
num_elements = numel(E_captured);
sigma_noise  = noise_abs_sum / (num_elements * sqrt(2/pi));
noise = sigma_noise * randn(size(E_captured));
E_captured = E_captured + noise;

%% 4. 光场层析重建：3层滑动窗口 Van Cittert 迭代（论文 3.3.2.3 节）
% ======================================================================
% 严格按论文式 3-34 到 3-40 的纯求和形式：
%   正向模型: E = Σ f_j * h_j         （式 3-34，无 Δy）
%   迭代更新: f^{k+1} = max(0, f^k + β·(E - H·f^k))  （式 3-38/39）
%   初值:     f_0 = E                  （式 3-40）
%
% 【问题12b】Van Cittert（论文式3-38原始形式，无H^T）
%   线性PSF → H 对称 → Van Cittert 收敛有保证
%   β = 1/N_win：DC模式一步精确收敛（|1-β·λ_DC|=|1-1|=0）
%   3层窗口: β=1/3, 边界2层窗口: β=1/2
%
% 关键区别：dy_sim 仅用于第3节的 400 层成像积分。
% 重建用 3 层，采用纯求和，无需 dy_rec。
% 纯求和体系下 f 的物理含义为"该层对图像的总辐射贡献"，
% 与辐射密度的关系为 f_sum = f_density · Δy_rec。
% 温度反解时需先除以 Δy_rec 还原为密度。
% ======================================================================
fprintf('\n开始3层滑动窗口 Van Cittert 重建 (20次迭代, β=1/N_win)...\n');

max_iter = 20;

f_est = zeros(Nz, Nx, N_rec);

% 中心层(y=0)温度真值
center_true_T = t1 * exp( - ( m_param * ( (Z_grid/Z_max).^2 + (X.^2)/(R^2) ) - n_param ).^2 ) + t2;

% --- 辐射→温度转换参数 ---
L_center = 25 / 1000;
rad_coeff = c1 * lambda^-5 * epsilon / (4*pi) * exp(-kappa * L_center);
rad2temp_raw = @(I_norm) c2 ./ (lambda * log(1 + rad_coeff ./ (max(I_norm, 1e-20) * I_max)));

% 辐射密度阈值（低于此值的像素视为背景）
f_thresh = 0.01;

% ======================================================================
% 温度域 SSIM：火焰 ROI 内计算
% ======================================================================
T_true_norm = (center_true_T - t2) / t1;

% 计算火焰 ROI（外接矩形）
fire_mask = center_true_T > (t2 + 50);  % T > 950K
[row_fire, col_fire] = find(fire_mask);
roi_r1 = min(row_fire); roi_r2 = max(row_fire);
roi_c1 = min(col_fire); roi_c2 = max(col_fire);
T_true_roi = T_true_norm(roi_r1:roi_r2, roi_c1:roi_c2);

fprintf('  火焰 ROI: 行[%d:%d], 列[%d:%d], 大小 %dx%d (全图 %dx%d)\n', ...
    roi_r1, roi_r2, roi_c1, roi_c2, ...
    roi_r2-roi_r1+1, roi_c2-roi_c1+1, Nz, Nx);
fprintf('  火焰像素占比: %.1f%%\n', 100*sum(fire_mask(:))/(Nz*Nx));

% 记录中心层 SSIM 曲线
ssim_roi_history  = zeros(max_iter, 1);
ssim_full_history = zeros(max_iter, 1);

% --- 逐层滑动窗口重建 ---
for target = 1:N_rec
    if target == 1
        win_idx = [1, 2];           target_pos_in_win = 1;
    elseif target == N_rec
        win_idx = [N_rec-1, N_rec]; target_pos_in_win = 2;
    else
        win_idx = [target-1, target, target+1]; target_pos_in_win = 2;
    end
    N_win = length(win_idx);

    % 【问题12b】β = 1/N_win（DC模式一步精确收敛）
    beta = 1 / N_win;

    fprintf('  目标层 %d (y=%.1fmm): 窗口层 [%s], 窗口大小 %d, β=%.3f\n', ...
        target, y_rec(target), num2str(win_idx), N_win, beta);

    % 4.1 构建窗口内 PSF 矩阵（线性插值σ，仅用于正向投影）
    PSF_win = cell(N_win, N_win);
    for ki = 1:N_win
        for ji = 1:N_win
            dy_abs = abs(y_rec(win_idx(ki)) - y_rec(win_idx(ji)));
            sigma = calc_sigma(dy_abs);
            psf_size = ceil(sigma * 3) * 2 + 1;
            PSF_win{ki, ji} = fspecial('gaussian', [psf_size, psf_size], sigma);
        end
    end

    % 4.2 初值 f_0 = E（论文式 3-40）
    f_win = zeros(Nz, Nx, N_win);
    for ji = 1:N_win
        f_win(:,:,ji) = E_captured(:,:,win_idx(ji));
    end

    % --- 初值诊断（仅中心层） ---
    if target == 4
        f0_center = f_win(:,:,target_pos_in_win);
        center_sim_idx = find(abs(y_sim) <= dy_rec/2);
        f_true_sum = sum(f_sim(:,:,center_sim_idx), 3) * dy_sim;
        fprintf('    [初值诊断] f_0: max=%.4f, f_true: max=%.4f, 比值=%.2fx\n', ...
            max(f0_center(:)), max(f_true_sum(:)), max(f0_center(:))/max(f_true_sum(:)));
    end

    % 4.3 Van Cittert 迭代（纯求和，无 dy_rec，无 H^T）
    for iter = 1:max_iter
        % 正向: E_est = Σ f_j * h_j（式 3-34，纯求和，无 Δy）
        E_est_win = zeros(Nz, Nx, N_win);
        for ki = 1:N_win
            for ji = 1:N_win
                E_est_win(:,:,ki) = E_est_win(:,:,ki) + ...
                    imfilter(f_win(:,:,ji), PSF_win{ki,ji}, 'replicate', 'same');
            end
        end

        % 计算残差 r_j = E_j - (Hf)_j
        residual_win = zeros(Nz, Nx, N_win);
        for ki = 1:N_win
            residual_win(:,:,ki) = E_captured(:,:,win_idx(ki)) - E_est_win(:,:,ki);
        end
        res_norm = norm(residual_win(:));

        % 计算更新量 Δf = β·r
        update_win = beta * residual_win;
        upd_norm = norm(update_win(:));

        % Van Cittert 更新: f_j += β·(E_j - (Hf)_j), 非负约束（式 3-38/39）
        f_win_old = f_win;
        for ji = 1:N_win
            f_win(:,:,ji) = max(0, f_win(:,:,ji) + update_win(:,:,ji));
        end

        % 实际变化量（含非负裁剪后）和裁剪比例
        actual_change = norm(f_win(:) - f_win_old(:));
        clipped = sum(f_win_old(:) + update_win(:) < 0) / numel(f_win);

        % 中心层 SSIM（火焰 ROI 内，温度域）+ 残差诊断
        if target == 4
            % 纯求和 f → 辐射密度（÷dy_rec）→ 温度（含阈值）→ 归一化 → SSIM
            f_density = f_win(:,:,target_pos_in_win) / dy_rec;
            T_iter = ones(Nz, Nx) * t2;
            valid = f_density > f_thresh;
            T_iter(valid) = rad2temp_raw(f_density(valid));
            T_iter = max(t2, min(t1 + t2, T_iter));
            T_iter_norm = (T_iter - t2) / t1;

            % ROI 内 SSIM
            T_iter_roi = T_iter_norm(roi_r1:roi_r2, roi_c1:roi_c2);
            iter_ssim = ssim(T_iter_roi, T_true_roi);

            % 全图 SSIM 作对照
            iter_ssim_full = ssim(T_iter_norm, T_true_norm);

            ssim_roi_history(iter)  = iter_ssim;
            ssim_full_history(iter) = iter_ssim_full;

            % E_est 与 E_captured 峰值比
            est_max = max(E_est_win(:));
            cap_max_in_win = 0;
            for ki = 1:N_win
                cap_max_in_win = max(cap_max_in_win, max(max(E_captured(:,:,win_idx(ki)))));
            end

            f_target = f_win(:,:,target_pos_in_win);
            fprintf('    Iter %2d/%d | ||r||=%.2f | ||βr||=%.4f | Δf=%.4f | clip=%.1f%% | f:[%.3f,%.3f] | E_est/E_cap=%.3f | SSIM: ROI=%.4f 全图=%.4f\n', ...
                iter, max_iter, res_norm, upd_norm, actual_change, ...
                100*clipped, min(f_target(:)), max(f_target(:)), ...
                est_max/cap_max_in_win, iter_ssim, iter_ssim_full);
        end
    end

    f_est(:,:,target) = f_win(:,:,target_pos_in_win);
end

% 最终评价
f_density_final = f_est(:,:,4) / dy_rec;
T_final = ones(Nz, Nx) * t2;
valid = f_density_final > f_thresh;
T_final(valid) = rad2temp_raw(f_density_final(valid));
T_final = max(t2, min(t1 + t2, T_final));
T_final_norm = (T_final - t2) / t1;
T_final_roi = T_final_norm(roi_r1:roi_r2, roi_c1:roi_c2);
ssim_final_roi  = ssim(T_final_roi, T_true_roi);
ssim_final_full = ssim(T_final_norm, T_true_norm);
fprintf('\n重建完成: 温度 SSIM: ROI=%.4f, 全图=%.4f\n', ssim_final_roi, ssim_final_full);

%% 5. 绘图：显示温度场（对应论文图4-5与图4-7）
figure('Name', '三维温度场重建结果 (400层模拟 -> 7层重建)', 'Position', [100, 100, 1200, 900]);

center_idx = 4;
T_true = center_true_T;

% 重建温度（含阈值处理）
f_density_recon = f_est(:,:,center_idx) / dy_rec;
T_recon = ones(Nz, Nx) * t2;
v1 = f_density_recon > f_thresh;
T_recon(v1) = rad2temp_raw(f_density_recon(v1));
T_recon = max(t2, min(t1 + t2, T_recon));

% 重聚焦原图温度（含阈值处理）
f_cap_density = E_captured(:,:,center_idx) / (N_sim * dy_sim);
T_captured = ones(Nz, Nx) * t2;
v2 = f_cap_density > f_thresh;
T_captured(v2) = rad2temp_raw(f_cap_density(v2));
T_captured = max(t2, min(t1 + t2, T_captured));

T_range = [900, 2200];

subplot(2, 3, 1);
imagesc(X_vec, Z_vec, T_true);
colormap(hot); colorbar; caxis(T_range);
set(gca, 'YDir', 'normal'); title('中心断层 温度真值 (y=0)');
xlabel('径向位置 X (mm)'); ylabel('轴向高度 Z (mm)'); axis image;

subplot(2, 3, 2);
imagesc(X_vec, Z_vec, T_captured);
colormap(hot); colorbar; caxis(T_range);
set(gca, 'YDir', 'normal'); title('重聚焦原图对应温度 (含层叠模糊)');
xlabel('径向位置 X (mm)'); axis image;

subplot(2, 3, 3);
imagesc(X_vec, Z_vec, T_recon);
colormap(hot); colorbar; caxis(T_range);
set(gca, 'YDir', 'normal'); title('Van Cittert 重建温度场 (20次迭代)');
xlabel('径向位置 X (mm)'); axis image;

% 画 ROI 框（绿色虚线标注火焰评价区域）
for sp = [1 2 3]
    subplot(2, 3, sp); hold on;
    roi_x = [X_vec(roi_c1), X_vec(roi_c2), X_vec(roi_c2), X_vec(roi_c1), X_vec(roi_c1)];
    roi_z = [Z_vec(roi_r1), Z_vec(roi_r1), Z_vec(roi_r2), Z_vec(roi_r2), Z_vec(roi_r1)];
    plot(roi_x, roi_z, 'g--', 'LineWidth', 1.5);
end

% 径向温度剖面对比
subplot(2, 3, 4);
z_target = 30;
[~, z_idx] = min(abs(Z_vec - z_target));
plot(X_vec, T_true(z_idx, :), 'k-', 'LineWidth', 2); hold on;
plot(X_vec, T_recon(z_idx, :), 'r--', 'LineWidth', 2);
plot(X_vec, T_captured(z_idx, :), 'b:', 'LineWidth', 1.5);
legend('温度真值', '重建温度', '含模糊原图温度', 'Location', 'best');
title(sprintf('高度 Z = %d mm 处径向温度剖面对比', z_target));
xlabel('径向位置 X (mm)'); ylabel('温度 (K)'); grid on;
ylim(T_range);

% SSIM 随迭代次数变化曲线（对应论文图 4-8(c)）
subplot(2, 3, 5);
plot(1:max_iter, ssim_roi_history, 'g-o', 'LineWidth', 2, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'g'); hold on;
plot(1:max_iter, ssim_full_history, 'b--s', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('迭代次数'); ylabel('SSIM');
title(sprintf('\\sigma_{neighbor} = %.1f, \\beta = 1/N_{win}', sigma_neighbor));
legend('ROI SSIM', '全图 SSIM', 'Location', 'southeast');
grid on; ylim([0 1]); xlim([1 max_iter]);

% 温度误差图
subplot(2, 3, 6);
T_error = T_recon - T_true;
imagesc(X_vec, Z_vec, T_error);
colorbar; caxis([-200, 200]);
set(gca, 'YDir', 'normal');
title('温度误差 (重建 - 真值) [K]');
xlabel('径向位置 X (mm)'); ylabel('轴向高度 Z (mm)'); axis image;

%% 6. 三维温度场可视化（m_vis=2.5 用于复现论文图4-5视觉效果）

figure('Name', '三维温度场可视化', 'Position', [150, 150, 1100, 500], 'Color', 'w');

m_vis = 2.5;

subplot(1, 2, 1);
hold on; grid on; box on;
[xc, yc, zc] = cylinder(R, 30);
zc = zc * Z_max;
mesh(xc, yc, zc, 'EdgeColor', [0.8 0 0], 'FaceColor', 'none', 'LineWidth', 0.5, 'FaceAlpha', 0);

T_center = t1 * exp( - ( m_vis * ( (Z_grid/Z_max).^2 + (X.^2 + 0^2)/(R^2) ) - n_param ).^2 ) + t2;
s_center = surf(X, zeros(size(X)), Z_grid, T_center, 'EdgeColor', 'none');
shading interp;
alpha_center = (T_center - 1000) / 100;
alpha_center(alpha_center > 1) = 1; alpha_center(alpha_center < 0) = 0;
set(s_center, 'FaceAlpha', 'interp', 'AlphaData', alpha_center);

view(35, 25);
colormap(hot); caxis([1000, 2200]);
xlabel('x'); ylabel('z'); zlabel('y');
title('(a)', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
set(gca, 'Color', [0.95 0.95 0.95]);
axis equal; xlim([-30 30]); ylim([-30 30]); zlim([0 Z_max]);
quiver3(0,0,0, 45,0,0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
quiver3(0,0,0, 0,45,0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
quiver3(0,0,0, 0,0,85, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);

subplot(1, 2, 2);
hold on; grid off; box on;

y_show = [-21, -14, -7, 0, 7, 14, 21];
for kk = 1:length(y_show)
    y_val = y_show(kk);
    T_slice = t1 * exp( - ( m_vis * ( (Z_grid/Z_max).^2 + (X.^2 + y_val^2)/(R^2) ) - n_param ).^2 ) + t2;
    s_slice = surf(X, y_val * ones(size(X)), Z_grid, T_slice, 'EdgeColor', 'none');
    shading interp;
    alpha_slice = (T_slice - 1000) / 100;
    alpha_slice(alpha_slice > 1) = 1; alpha_slice(alpha_slice < 0) = 0;
    set(s_slice, 'FaceAlpha', 'interp', 'AlphaData', alpha_slice);
    plot3([-21 21 21 -21 -21], [y_val y_val y_val y_val y_val], ...
          [0 0 Z_max Z_max 0], 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    text(0, y_val, Z_max + 4, sprintf('%.0f', y_val), ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

view(-45, 25);
colormap(hot); caxis([1000, 2200]);
xlabel('x'); ylabel('Depth (z)'); zlabel('y');
title('(b)', 'Units', 'normalized', 'Position', [0.5, -0.15, 0]);
xlim([-21 21]); ylim([-25 25]); zlim([0 Z_max]);
pbaspect([1 2.5 1.2]);
set(gca, 'ZTick', [], 'XTick', [], 'YDir', 'reverse');

cb = colorbar('southoutside');
cb.Position = [0.65 0.15 0.25 0.03];
cb.Ticks = 1000:200:2200;
title(cb, 'Temperature [^\circC]', 'FontSize', 11);