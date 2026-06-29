function burstTrain = DetectBurstsNeuro(SpikeTimesInSamples, SampFreq, RecDuration, Params)
% DetectBurstsNeuro  Single-channel burst detection via the logarithmic
%   ISI-histogram method (construct_ISIhist + calcISImax + burstDetection).
%
% INPUT
%   SpikeTimesInSamples : spike times of one channel, in samples
%   SampFreq            : sampling frequency [Hz]
%   RecDuration         : recording duration [s]
%   Params              : struct with field Log_ISI_Hist containing:
%                           BinsperDec, VoidThr, DefaultThr, MaxThr,
%                           MinSpikes, Smoother, Span, MPD
%
% OUTPUT
%   burstTrain : [nBursts x 6+] burst matrix from burstDetection (empty if
%                fewer than N spikes). See burstDetection for column layout.

sf          = SampFreq;
nBinsPerDec = Params.Log_ISI_Hist.BinsperDec;
voidThr     = Params.Log_ISI_Hist.VoidThr;
ISIThr      = Params.Log_ISI_Hist.DefaultThr;
maxISIThr   = Params.Log_ISI_Hist.MaxThr;
minNumSpks  = Params.Log_ISI_Hist.MinSpikes;
smoother    = Params.Log_ISI_Hist.Smoother;
span        = Params.Log_ISI_Hist.Span;
mpd         = Params.Log_ISI_Hist.MPD;
N           = 2;   % original single-channel algorithm uses N = 2

if numel(SpikeTimesInSamples) > N
    totaldur = RecDuration*SampFreq;
    peakTrn  = zeros(totaldur, 1);
    peakTrn(SpikeTimesInSamples) = 1;
    sample = SpikeTimesInSamples/sf;

    [bins, ~, histogram] = construct_ISIhist(sample, N, span, nBinsPerDec, smoother);
    [ISImax, ~, flags]   = calcISImax(bins, histogram, mpd, voidThr, ISIThr);
    [burstTrain, ~]      = burstDetection(peakTrn, ISImax, flags, ISIThr, maxISIThr, sf, minNumSpks);
else
    burstTrain = [];
end
end
