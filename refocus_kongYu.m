function refocused_Im = refocus_kongYu(LF, alpha)

[U, V, row, col, ~] = size(LF);
M = U * V;

LF_Image = zeros(M, row, col, 3, 'double');

% 把每个子孔径图像取出来
for u = 1:U
    for v = 1:V
        k = (u-1)*V + v;
        LF_Image(k, :, :, :) = double(squeeze(LF(u, v, :, :, :)));
    end
end

% 构造标准的视角坐标
u_coords = ((1:U) - (U+1)/2) * alpha;
v_coords = ((1:V) - (V+1)/2) * alpha;

% 注意这里用 ndgrid，保证和 (u,v) 顺序一致
[UU, VV] = ndgrid(u_coords, v_coords);

D = zeros(M, 2);
for u = 1:U
    for v = 1:V
        k = (u-1)*V + v;
        D(k,:) = [UU(u,v), VV(u,v)];
    end
end

shifted_Imr = zeros(M, row, col);
shifted_Img = zeros(M, row, col);
shifted_Imb = zeros(M, row, col);

for k = 1:M
    shifted_Imr(k,:,:) = ImWarp(squeeze(LF_Image(k,:,:,1)), -D(k,1), -D(k,2));
    shifted_Img(k,:,:) = ImWarp(squeeze(LF_Image(k,:,:,2)), -D(k,1), -D(k,2));
    shifted_Imb(k,:,:) = ImWarp(squeeze(LF_Image(k,:,:,3)), -D(k,1), -D(k,2));
end

refocused_Im = zeros(row, col, 3);
refocused_Im(:,:,1) = squeeze(mean(shifted_Imr, 1));
refocused_Im(:,:,2) = squeeze(mean(shifted_Img, 1));
refocused_Im(:,:,3) = squeeze(mean(shifted_Imb, 1));

end