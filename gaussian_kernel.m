function kernel = gaussian_kernel(kernel_size, sigma)
% gaussian_kernel  Build a 1-D Gaussian smoothing kernel.
%
% INPUT
%   kernel_size : kernel width in samples
%   sigma       : standard deviation in samples
%
% OUTPUT
%   kernel      : column vector with the Gaussian kernel (not normalised)

x = (-(kernel_size-1)/2 : (kernel_size-1)/2)';
kernel = exp(-x.^2 / (2*sigma^2)) / (sigma * sqrt(2*pi));
end
