function [NBs, peaks, smoothed] = detectNetworkBursts_Doorn2024(spikeratetot, APsinbin, timeBin, movingWindow_sec, sigma, kernelsize)
% detectNetworkBursts_Doorn2024  Network-burst detection on a smoothed global
%   spike-rate signal with dynamic, locally-adaptive start/stop thresholds
%   (van der Molen / Doorn 2024 approach).
%
% INPUT
%   spikeratetot     : global spike-rate signal (spikes per bin), 1-by-nBins
%   APsinbin         : binned spike counts per electrode, [nElectrodes x nBins]
%   timeBin          : bin width [s]
%   movingWindow_sec : moving window for the local maximum [s]
%   sigma            : Gaussian smoothing std [samples]
%   kernelsize       : Gaussian kernel width [samples]
%
% OUTPUT
%   NBs      : [nNB x 5] matrix:
%                1 onset [bin], 2 end [bin], 3 total spikes in burst,
%                4 unused, 5 number of peaks inside the burst
%   peaks    : peak locations on the smoothed signal [bin]
%   smoothed : the smoothed global spike-rate signal

kernel   = gaussian_kernel(kernelsize, sigma);
kernel   = kernel / sum(kernel);
smoothed = smooth_signal(spikeratetot, kernel);

movingint = round(movingWindow_sec / timeBin);   % moving window [bins]
Tth       = round(50e-3 / timeBin);              % 50 ms [bins]

actelec       = sum(mean(APsinbin,2) > 0.02);    % number of active electrodes
Eth           = 0.3 * actelec;                   % min electrode participation
MBth          = max(spikeratetot) / 16;          % min peak height
prominence_th = max(spikeratetot) / 10;          % min peak prominence

[~, peaks] = findpeaks(smoothed, ...
    'MinPeakHeight', MBth, ...
    'MinPeakProminence', prominence_th);

NBs     = zeros(1000,5);
NBcount = 0;
i       = 1;
spikeratetot = smoothed;   % run the detector on the smoothed signal

while (i + Tth) <= length(spikeratetot)

    % Local maximum around the current bin -> dynamic thresholds
    if i > movingint && i < length(spikeratetot) - movingint
        localmax = max(spikeratetot(i-movingint : i+movingint));
    elseif i < length(spikeratetot) - movingint
        localmax = max(spikeratetot(1 : i+movingint));
    else
        localmax = max(spikeratetot(i-movingint : end));
    end

    Startth = 0.25 * localmax;
    Stopth  = 0.02 * localmax;

    % Burst start: signal above Startth for the whole 50 ms window AND enough
    % electrodes active in that same window.
    if all(spikeratetot(i:i+Tth-1) > Startth) && ...
       sum(sum(APsinbin(:,i:i+Tth-1),2) > Tth) > Eth

        NBcount = NBcount + 1;
        NBs(NBcount,1) = i;
        NBs(NBcount,3) = sum(sum(APsinbin(:,i:i+Tth-1)));
        i = i + Tth;

        % Extend the burst while activity stays above the stop criterion
        while i + 2*Tth - 1 <= length(spikeratetot) && ...
              any(spikeratetot(i:i+2*Tth-1) > Stopth)
            NBs(NBcount,3) = NBs(NBcount,3) + sum(APsinbin(:,i));
            i = i + 1;
        end

        NBs(NBcount,2) = i - 1;
        NBs(NBcount,5) = sum((peaks > NBs(NBcount,1)) & (peaks <= NBs(NBcount,2)));
    else
        i = i + 1;
    end
end

NBs = NBs(1:NBcount,:);
end
