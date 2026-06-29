function NB = DetectNetworkBurstsNeuroV2(app, Spikes, Params)
% DetectNetworkBurstsNeuroV2  Network-burst detection (Tampere / log ISI-N),
%   V2b threshold logic: peaks are selected with findpeaks using both
%   MinPeakHeight and MinPeakProminence (10% of max), and the ISI-N threshold
%   is the global minimum of the valley between the first two peaks. This
%   avoids spurious smoothing ripples being chosen as the second peak.
%   Merge / cut / channel filters / fragment analysis match DetectNetworkBurstsNeuro.
%
% INPUT
%   app    : app/object handle ([] -> mock struct); optional metadata fields
%            cellLineName, meaPlate, DIV, wellName
%   Spikes : table for one well with columns Time [s] and Channel
%   Params : struct with fields
%             .LogISINHist.{Smoother,BinsPerDec,Span,N}
%             .MinCount, .MinDuration, .MergingType, .MergingFracOrThres,
%             .MinChannels, .MinDomChannels, .MinChannelContribution
%             .timeBin (def 25e-3), .kernelsize (def 13), .sigma (def 4)
%
% OUTPUT
%   NB : struct with accepted NB .starts/.ends [s], per-burst channel tables
%        (.CHANNELSFIN/.DOMCHANNFIN), .ISINThreshold [ms], fragment fields
%        (.num_fragments, .isfragmented, .fragmentPeakTimes, .peakTimes) and
%        activity traces. Empty structure (WriteNoNB) when no valid NB found.

if isempty(app)
    app = struct();
end
if ~isfield(app,'cellLineName'), app.cellLineName = '?'; end
if ~isfield(app,'meaPlate'),     app.meaPlate     = '?'; end
if ~isfield(app,'DIV'),          app.DIV          = '?'; end
if ~isfield(app,'wellName'),     app.wellName     = '?'; end

NB = WriteNoNB();

% ISI-N histogram
Spikes     = table2struct(Spikes, 'ToScalar', true);
SpikeTimes = Spikes.Time;

[Steps, ~, HH] = construct_ISIhist(SpikeTimes, ...
    Params.LogISINHist.N, Params.LogISINHist.Span, ...
    Params.LogISINHist.BinsPerDec, Params.LogISINHist.Smoother);

if isempty(Steps) || isempty(HH)
    return
end

HHrnd = round(HH, 9);

% Threshold = global minimum of the valley between the first two relevant peaks
minPeakH      = 0.10 * max(HHrnd);
minProminence = 0.10 * max(HHrnd);
[~, peakIdx] = findpeaks(double(HHrnd), ...
    'MinPeakHeight', minPeakH, 'MinPeakProminence', minProminence);

if numel(peakIdx) < 2
    return
end

p1 = peakIdx(1);
p2 = peakIdx(2);
if p2 <= p1 + 1
    return
end

segment      = HHrnd(p1+1 : p2-1);
[~, relIdx]  = min(segment);
ThresholdInd = p1 + relIdx;

ISINThreshold_ms = Steps(ThresholdInd);
ISINThreshold_s  = ISINThreshold_ms / 1000;

[Burst, ~, ~] = BurstDetectISIn_MOD_AV(Spikes, Params.LogISINHist.N, ISINThreshold_s);

% General NB criteria
if numel(Burst.T_start) < Params.MinCount || isempty(Burst.T_start)
    return
end

Merging_Parameter = Params.MergingFracOrThres;

if numel(Burst.T_start) > 1 && strcmp(Params.MergingType, 'Fixed')
    [nbstartsmerged, nbendsmerged, CHANNELS_merged] = ...
        MergeBurstsUniversal(app, Burst.T_start, Burst.T_end, Merging_Parameter, 'Static', Burst.CN);
elseif numel(Burst.T_end) > 1 && strcmp(Params.MergingType, 'Adaptive')
    [nbstartsmerged, nbendsmerged, CHANNELS_merged] = ...
        MergeBurstsUniversal(app, Burst.T_start, Burst.T_end, Merging_Parameter, 'Adaptive', Burst.CN);
elseif numel(Burst.T_end) == 1 || strcmp(Params.MergingType, 'None')
    nbendsmerged    = Burst.T_end;
    nbstartsmerged  = Burst.T_start;
    CHANNELS_merged = Burst.CN;
else
    nbendsmerged    = [];
    nbstartsmerged  = [];
    CHANNELS_merged = {};
end

if isempty(nbendsmerged) || numel(nbendsmerged) < Params.MinCount
    return
end

% MinDuration cut
if Params.MinDuration > 0
    [nbstartscut, nbendscut, CHANNELS_cut] = ...
        CutBurstsUniversal(app, nbstartsmerged, nbendsmerged, Params.MinDuration, CHANNELS_merged);
else
    nbstartscut  = nbstartsmerged;
    nbendscut    = nbendsmerged;
    CHANNELS_cut = CHANNELS_merged;
end

if isempty(nbendscut) || numel(nbendscut) < Params.MinCount
    return
end

% MinChannels filter
NCH = cell(1, numel(CHANNELS_cut));
for k = 1:numel(CHANNELS_cut)
    NCH{1,k} = unique(CHANNELS_cut{k});
end
BurstChannels        = cellfun(@numel, NCH(:))';
worthy_burst_indices = find(BurstChannels(:) >= Params.MinChannels);

nbstarts_MC = nbstartscut(worthy_burst_indices);
nbends_MC   = nbendscut(worthy_burst_indices);
CHANNELS_MC = CHANNELS_cut(1, worthy_burst_indices);

if isempty(nbends_MC) || numel(nbends_MC) < Params.MinCount
    return
end

% MinChannelContribution filter
CHANNELS_table    = cell(1, length(CHANNELS_MC));
mean_contribution = zeros(1, length(CHANNELS_MC));
for bnmb = 1:length(CHANNELS_MC)
    Electrode = unique(CHANNELS_MC{bnmb});
    mean_contribution(bnmb) = numel(CHANNELS_MC{bnmb}) / numel(Electrode);
    SpikeCount = zeros(numel(Electrode),1);
    for chnl = 1:numel(Electrode)
        SpikeCount(chnl,1) = nnz(strcmp(CHANNELS_MC{bnmb}, Electrode{chnl}));
    end
    CHANNELS_table{bnmb} = table(Electrode, SpikeCount);
    clear Electrode SpikeCount
end

CHANNELS4 = CHANNELS_table;
for iiii = 1:length(CHANNELS4)
    toDelete = CHANNELS4{iiii}.SpikeCount < ...
        floor(Params.MinChannelContribution * mean_contribution(iiii));
    CHANNELS4{iiii}(toDelete,:) = [];
    clear toDelete
end

contrchannelsthr      = cellfun(@height, CHANNELS4);
worthy_burst_indices2 = find(contrchannelsthr >= Params.MinDomChannels);

if numel(worthy_burst_indices2) < Params.MinCount
    return
end

NB.starts        = round(nbstarts_MC(worthy_burst_indices2), 6);
NB.ends          = round(nbends_MC(worthy_burst_indices2), 6);
NB.CHANNELSFIN   = CHANNELS_table(1, worthy_burst_indices2);
NB.DOMCHANNFIN   = CHANNELS4(1, worthy_burst_indices2);
NB.ISINThreshold = ISINThreshold_ms;

% Fragment (sub-burst) analysis
if ~isfield(Params,'timeBin'),    Params.timeBin    = 25e-3; end
if ~isfield(Params,'kernelsize'), Params.kernelsize = 13;    end
if ~isfield(Params,'sigma'),      Params.sigma      = 4;     end

timeBin   = Params.timeBin;
startTime = 0;
stopTime  = 600;
rectime   = stopTime - startTime;

Channel = Spikes.Channel;
Time    = Spikes.Time;

validMask   = (Time > startTime) & (Time < stopTime);
Time_wot    = Time(validMask) - startTime;
Channel_wot = Channel(validMask);

electrodes  = unique(Channel_wot);
nElectrodes = numel(electrodes);
nBins       = floor(rectime / timeBin);

APsinbin  = zeros(nElectrodes, nBins);
spikerate = zeros(nElectrodes, nBins);
for e = 1:nElectrodes
    APt = Time_wot(strcmp(Channel_wot, electrodes{e}));
    for k = 1:nBins
        t0 = (k-1) * timeBin;
        t1 =  k    * timeBin;
        APsinbin(e,k)  = sum((APt > t0) & (APt < t1));
        spikerate(e,k) = APsinbin(e,k) / timeBin;
    end
end

APsinbintot           = sum(APsinbin,  1);
spikeratetot          = sum(spikerate, 1);
activeChannelsPerBin  = sum(APsinbin > 0, 1);
spikeratetot_weighted = spikeratetot .* activeChannelsPerBin;

windowBins_plot = ceil(0.1 / timeBin);
spikerateplot   = movmean(spikeratetot, windowBins_plot);

kernel = gaussian_kernel(Params.kernelsize, Params.sigma);
kernel = kernel / sum(kernel);
spikeratesmooth = smooth_signal(spikeratetot, kernel);

MBth          = max(spikeratetot) / 16;
prominence_th = max(spikeratetot) / 10;
[~, locs] = findpeaks(spikeratesmooth, ...
    'MinPeakHeight', MBth, 'MinPeakProminence', prominence_th);

peakTimes = (locs - 1) * timeBin;

NB.num_fragments     = zeros(size(NB.starts));
NB.fragmentPeakTimes = cell(size(NB.starts));
for i = 1:length(NB.starts)
    fragment_mask           = (peakTimes >= NB.starts(i)) & (peakTimes < NB.ends(i));
    NB.num_fragments(i)     = sum(fragment_mask);
    NB.fragmentPeakTimes{i} = peakTimes(fragment_mask);
end

NB.isfragmented            = NB.num_fragments >= 2;
NB.peakTimes               = peakTimes;
NB.spikeTimes              = (0:length(spikeratesmooth)-1) * timeBin;
NB.spikeratetot_frag       = spikeratetot;
NB.spikeratesmooth         = spikeratesmooth;
NB.spikerateplot           = spikerateplot;
NB.APsinbinfrag            = APsinbin;
NB.APsinbintotfrag         = APsinbintot;
NB.activeChannelsPerBin    = activeChannelsPerBin;
NB.spikeratetotweighted    = spikeratetot_weighted;
NB.spikerateplotweighted   = movmean(spikeratetot_weighted, windowBins_plot);
NB.spikeratesmoothweighted = smooth_signal(spikeratetot_weighted, kernel);
end


function NB = WriteNoNB()
% Default empty output used when no valid network burst is found.
NB.starts        = [];
NB.ends          = [];
NB.CHANNELSFIN   = [];
NB.DOMCHANNFIN   = [];
NB.ISINThreshold = nan;
end
