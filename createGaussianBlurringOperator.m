function [Op, Op_adj, fftFilter] = createGaussianBlurringOperator(size_I,hsize,sigma)

% Creates a operator for a blurring with a Gaussian kernel of size 
% 'hsize' and standard deviation 'sigma'

kernel = fspecial('gaussian', hsize, sigma);

fftFilter = zeros(size_I);
fftFilter(1:size(kernel,1),1:size(kernel,2)) = kernel;

% Center the kernel
fftFilter = circshift(fftFilter,-(size(kernel,1)-1)/2,1);
fftFilter = circshift(fftFilter,-(size(kernel,2)-1)/2,2);

% Precalculate FFT
fftFilter = fftn(fftFilter);
fftFilterC = conj(fftFilter);

% Setup the operators
Op = @(x) ifftn(fftn(reshape(x,size_I)) .* fftFilter);
Op_adj = @(y) ifftn(fftn(reshape(y,size_I)) .* fftFilterC);

end