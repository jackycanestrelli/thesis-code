function spikeratesmooth = smooth_signal(x, kernel)
% smooth_signal  Convolve a signal with a kernel (same-length output).
%
% INPUT
%   x      : input signal (vector)
%   kernel : smoothing kernel (vector), e.g. from gaussian_kernel
%
% OUTPUT
%   spikeratesmooth : smoothed signal, same length as x

spikeratesmooth = conv(x, kernel, 'same');
end
