function recon = deconvRL_operator(forwardFUN, backwardFUN, img, num_iter, init)

recon = init;

fprintf('\nDeconvolution:\n');

for i = 1:num_iter
    tic;

    fpj = forwardFUN(recon);
    fpj = real(fpj);

    errorBack = img ./ (fpj + 1e-10);

    errorBack(isnan(errorBack)) = 0;
    errorBack(isinf(errorBack)) = 0;
    errorBack(errorBack < 0) = 0;
    errorBack(errorBack > 1e10) = 0;

    bpjError = backwardFUN(errorBack);
    bpjError = real(bpjError);

    recon = recon .* bpjError;
    recon(recon < 0) = 0;

    elapsed = toc;
    fprintf('iter %d | %d, took %.4f secs\n', i, num_iter, elapsed);
end

end