function recon = deconvRL_operator2(forwardFUN, backwardFUN, img, num_iter, init)

recon = init;

fprintf('\nDeconvolution:\n');

% RL 归一化项：A^T * 1
normalizer = backwardFUN(ones(size(img)));
normalizer = real(normalizer);
normalizer(normalizer <= 0) = 1e-10;

for i = 1:num_iter
    tic;

    % 正向投影：A * x_k
    fpj = forwardFUN(recon);
    fpj = real(fpj);

    % 比例误差：b ./ (A*x_k)
    errorBack = img ./ (fpj + 1e-10);

    % 数值稳定处理
    errorBack(isnan(errorBack)) = 0;
    errorBack(isinf(errorBack)) = 0;
    errorBack(errorBack < 0) = 0;
    errorBack(errorBack > 1e10) = 0;

    % 反向投影：A^T * errorBack
    bpjError = backwardFUN(errorBack);
    bpjError = real(bpjError);

    % RL 更新：
    % x_{k+1} = x_k .* A^T(b ./ A*x_k) ./ A^T(1)
    recon = recon .* bpjError ./ normalizer;

    % 非负约束
    recon(recon < 0) = 0;
    recon(isnan(recon)) = 0;
    recon(isinf(recon)) = 0;

    elapsed = toc;
    fprintf('iter %d | %d, took %.4f secs\n', i, num_iter, elapsed);
end

end