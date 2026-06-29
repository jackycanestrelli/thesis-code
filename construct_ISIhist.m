function [Steps, ISI_N, histogram] = construct_ISIhist(SpikeTimes, N, span, nBinsPerDec, smoother)
% construct_ISIhist  Build a smoothed, normalised logarithmic ISI_N histogram.
%
% INPUT
%   SpikeTimes  : spike times of one channel [s], ascending
%   N           : ISI order (distance between the i-th and (i+N-1)-th spike)
%   span        : smoothing span passed to smooth()
%   nBinsPerDec : number of histogram bins per decade
%   smoother    : smoothing method name passed to smooth() (e.g. 'moving')
%
% OUTPUT
%   Steps     : logarithmic bin edges of the ISI_N histogram [ms]
%   ISI_N     : raw ISI_N values [ms]
%   histogram : normalised + smoothed histogram (column vector, >= 0)
%
% Returns empty outputs when there are not more than N spikes.
% Note: bins start at 1 ms (logspace from 0); intended for N >= 10.

if numel(SpikeTimes) > N
    ISI_N = (SpikeTimes( N:end ) - SpikeTimes( 1:end-(N-1) ))*1000;
    Steps = logspace(0, ceil(max(log10(ISI_N))), nBinsPerDec*ceil(max(log10(ISI_N))));

    counts = histcounts(ISI_N, Steps);
    % Add an extra bin to mimic the legacy histc() behaviour (count of the right edge).
    counts = [counts sum(ISI_N==Steps(end))];

    norm_counts = counts/sum(counts);
    histogram = smooth(norm_counts, span, smoother);
    histogram(histogram<0) = 0;
    histogram = histogram(:);
else
    Steps = [];
    ISI_N = [];
    histogram = [];
end
end
