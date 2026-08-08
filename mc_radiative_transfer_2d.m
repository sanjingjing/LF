function I_img = mc_radiative_transfer_2d(Eb, hot_mask, y, z, p)
% =========================================================
% Monte Carlo 2D 火焰介质辐射传输
%
% 依据你的截图公式：
%   (1) 发射：由 Eb 给出
%   (2) 自由程：Li = -ln(R)/ke
%   (3) 碰撞后按 omega 判断吸收/散射
%   (4) 从观察侧（y = ymin）逃逸才记入图像
%
% 输入：
%   Eb        - 单色黑体辐射场
%   hot_mask  - 发射区域掩膜
%   y, z      - 物理网格
%   p.ke      - 消光系数
%   p.omega   - 散射反照率
%   p.N_rays  - 每个发射体元发射光线数
%   p.max_scatter
%   p.use_cone_emission
%   p.alpha_deg
%   p.use_anisotropic
%   p.A1
%   p.ymin/ymax/zmin/zmax
%
% 输出：
%   I_img     - 观测面有效辐射图
% =========================================================

I_img = zeros(size(Eb));
hot_idx = find(hot_mask);

alpha = deg2rad(p.alpha_deg);

for n = 1:length(hot_idx)
    idx = hot_idx(n);
    [r, c] = ind2sub(size(Eb), idx);

    E0 = Eb(r,c);
    if E0 <= 0
        continue;
    end

    y0 = y(r,c);
    z0 = z(r,c);

    % 每条射线携带的能量
    ray_energy = E0 / p.N_rays;

    escaped_energy = 0;

    for k = 1:p.N_rays
        pos = [y0, z0];

        % -------------------------
        % 初始方向采样
        % -------------------------
        if p.use_cone_emission
            % 朝向主镜头的小立体角发射
            % 这里把主镜头方向近似为 y 负方向
            theta = pi + (2*rand - 1) * alpha;
        else
            % 半空间发射
            theta = pi/2 + pi*rand;
        end

        dir = [cos(theta), sin(theta)];
        alive = true;
        scatter_count = 0;

        while alive
            % -------------------------
            % 自由程采样：Li = -ln(R)/ke
            % -------------------------
            Li = -log(max(rand,1e-12)) / p.ke;
            new_pos = pos + Li * dir;

            % -------------------------
            % 边界判断
            % -------------------------
            if new_pos(1) <= p.ymin
                % 从观测侧逃逸，计入有效辐射
                escaped_energy = escaped_energy + ray_energy;
                break;
            end

            if new_pos(1) >= p.ymax || new_pos(2) <= p.zmin || new_pos(2) >= p.zmax
                % 从其他边界逃逸，不计入镜头
                break;
            end

            % -------------------------
            % 碰撞后吸收/散射判定
            % -------------------------
            if rand > p.omega
                % 吸收
                alive = false;
            else
                % 散射
                scatter_count = scatter_count + 1;
                if scatter_count > p.max_scatter
                    alive = false;
                    break;
                end

                % 更新位置到碰撞点
                pos = new_pos;

                % -------------------------
                % 散射方向采样
                % -------------------------
                if p.use_anisotropic
                    theta_new = sample_theta_A1(p.A1);
                else
                    % 各向同性散射
                    theta_new = 2*pi*rand;
                end

                dir = [cos(theta_new), sin(theta_new)];
            end
        end
    end

    I_img(r,c) = escaped_energy;
end

% 为避免数值过大，按每体元射线数归一化
I_img = I_img / p.N_rays;

end

% =========================================================
function theta = sample_theta_A1(A1)
% 一阶相函数近似：
%   Phi(theta) = 1 + A1*cos(theta)
% 用拒绝采样实现
%
% A1 = 0  -> 各向同性
% A1 > 0  -> 前向散射偏强
% A1 < 0  -> 后向散射偏强

if abs(A1) < 1e-12
    theta = 2*pi*rand;
    return;
end

max_pdf = 1 + abs(A1);

while true
    theta_try = 2*pi*rand;
    pdf_val = 1 + A1*cos(theta_try);
    if pdf_val < 0
        pdf_val = 0;
    end

    if rand * max_pdf <= pdf_val
        theta = theta_try;
        return;
    end
end
end