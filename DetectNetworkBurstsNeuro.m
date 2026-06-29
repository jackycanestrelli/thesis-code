function NB = DetectNetworkBurstsNeuro(app, Spikes, Params)
% DetectNetworkBurstsNeuro  Network-burst detection (Tampere / log ISI-N method).
%   Builds the pooled log ISI-N histogram, picks the ISI-N threshold from the
%   first local minimum after the main peak, detects candidate bursts, then
%   applies merge / duration / channel-count / channel-contribution criteria
%   and a fragment (sub-burst) analysis on the accepted network bursts.
%
% INPUT
%   app    : app/object handle passed to Merge/Cut helpers ([] -> mock struct)
%   Spikes : table for one well with columns Time [s] and Channel
%   Params : struct with fields
%             .LogISINHist.{Smoother,BinsPerDec,Span,N}
%             .MinCount, .MinDuration, .MergingType ('Fixed'|'Adaptive'|'None'),
%             .MergingFracOrThres, .MinChannels, .MinDomChannels,
%             .MinChannelContribution, .timeBin, .kernelsize, .sigma
%
% OUTPUT
%   NB : struct with fields
%          .starts, .ends                accepted NB start/end times [s]
%          .CHANNELSFIN, .DOMCHANNFIN    per-burst channel/dominant-channel tables
%          .num_fragments, .is_fragmented, .fragmentPeakTimes, .peakTimes
%          activity traces (.spikeratetot_frag, .spikeratesmooth, ...)
%        Empty structure (via WriteNoNB) when no valid NB is found.

if isempty(app)
    app = struct();
end

NB = WriteNoNB();

% Convert spikes table into a scalar struct for downstream functions
Spikes = table2struct(Spikes, "ToScalar", true);
SpikeTimes = Spikes.Time;

% Build the log ISI-N histogram
[Steps, ~, HH] = construct_ISIhist(SpikeTimes, Params.LogISINHist.N, ...
    Params.LogISINHist.Span, Params.LogISINHist.BinsPerDec, ...
    Params.LogISINHist.Smoother);

% Threshold = first local minimum after the global maximum of the histogram
HHrnd = round(HH, 9);
[~, MaxInd] = max(HHrnd);
Localmins = find(islocalmin(HHrnd));
Localmins = Localmins(Localmins > MaxInd);

if ~isempty(Localmins)
    ThresholdInd     = Localmins(1);
    ISINThreshold_ms = Steps(ThresholdInd);
    ISINThreshold_s  = ISINThreshold_ms / 1000;
    [Burst, ~, ~]    = BurstDetectISIn_MOD_AV(Spikes, Params.LogISINHist.N, ISINThreshold_s);
else
    NB = WriteNoNB();
    return
end

% General network-burst criteria
if numel(Burst.T_start) < Params.MinCount || isempty(Burst.T_start)
    NB = WriteNoNB();
    return
end

Merging_Parameter = Params.MergingFracOrThres;

% Merge candidate bursts
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
end

if ~isempty(nbendsmerged) && numel(nbendsmerged) >= Params.MinCount

    % Minimum-duration cut
    if Params.MinDuration > 0
        [nbstartscut, nbendscut, CHANNELS_cut] = ...
            CutBurstsUniversal(app, nbstartsmerged, nbendsmerged, Params.MinDuration, CHANNELS_merged);
    else
        nbstartscut  = nbstartsmerged;
        nbendscut    = nbendsmerged;
        CHANNELS_cut = CHANNELS_merged;
    end

    if ~isempty(nbendscut) && numel(nbendscut) >= Params.MinCount

        % Minimum number of participating channels
        NCH = cell(1, numel(CHANNELS_cut));
        for k = 1:numel(CHANNELS_cut)
            NCH{1,k} = unique(CHANNELS_cut{k});
        end
        BurstChannels = cellfun(@numel, NCH(:))';
        worthy_burst_indices = find(BurstChannels(:) >= Params.MinChannels);

        nbstarts_MC = nbstartscut(worthy_burst_indices);
        nbends_MC   = nbendscut(worthy_burst_indices);
        CHANNELS_MC = CHANNELS_cut(1, worthy_burst_indices);

        % Minimum channel-contribution criterion
        if ~isempty(nbends_MC) && numel(nbends_MC) >= Params.MinCount

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

            % Keep electrodes above MinChannelContribution * mean contribution
            CHANNELS4 = CHANNELS_table;
            for iiii = 1:length(CHANNELS4)
                toDelete = CHANNELS4{iiii}.SpikeCount < ...
                    floor(Params.MinChannelContribution * mean_contribution(iiii));
                CHANNELS4{iiii}(toDelete,:) = [];
                clear toDelete
            end

            contrchannelsthr = cellfun(@height, CHANNELS4);
            worthy_burst_indices2 = find(contrchannelsthr >= Params.MinDomChannels);

            if numel(worthy_burst_indices2) >= Params.MinCount

                NB.starts      = round(nbstarts_MC(worthy_burst_indices2), 6);
                NB.ends        = round(nbends_MC(worthy_burst_indices2), 6);
                NB.CHANNELSFIN = CHANNELS_table(1, worthy_burst_indices2);
                NB.DOMCHANNFIN = CHANNELS4(1, worthy_burst_indices2);

                % --- Fragment (sub-burst) analysis ---
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
                        t1 = k * timeBin;
                        APsinbin(e,k)  = sum((APt > t0) & (APt < t1));
                        spikerate(e,k) = APsinbin(e,k) / timeBin;
                    end
                end

                APsinbintot          = sum(APsinbin, 1);
                spikeratetot         = sum(spikerate, 1);
                activeChannelsPerBin = sum(APsinbin > 0, 1);
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

                spikeTimesAxis = (0:length(spikeratesmooth)-1) * timeBin;
                peakTimes      = spikeTimesAxis(locs);

                % Count fragments (peaks) inside each accepted NB
                NB.num_fragments     = zeros(size(NB.starts));
                NB.fragmentPeakTimes = cell(size(NB.starts));
                for i = 1:length(NB.starts)
                    fragment_mask = (peakTimes >= NB.starts(i)) & (peakTimes < NB.ends(i));
                    NB.num_fragments(i)     = sum(fragment_mask);
                    NB.fragmentPeakTimes{i} = peakTimes(fragment_mask);
                end

                NB.is_fragmented = NB.num_fragments >= 2;
                NB.peakTimes     = peakTimes;
                NB.spikeTimes    = spikeTimesAxis;

                NB.spikeratetot_frag        = spikeratetot;
                NB.spikeratesmooth          = spikeratesmooth;
                NB.spikerateplot            = spikerateplot;
                NB.APsinbin_frag            = APsinbin;
                NB.APsinbintot_frag         = APsinbintot;
                NB.activeChannelsPerBin     = activeChannelsPerBin;
                NB.spikeratetot_weighted    = spikeratetot_weighted;
                NB.spikerateplot_weighted   = movmean(spikeratetot_weighted, windowBins_plot);
                NB.spikeratesmooth_weighted = smooth_signal(spikeratetot_weighted, kernel);
            else
                NB = WriteNoNB();
            end
        else
            NB = WriteNoNB();
        end
    else
        NB = WriteNoNB();
    end
else
    NB = WriteNoNB();
end
end


function NB = WriteNoNB()
% Default empty output used when no valid network burst is found.
NB.starts        = [];
NB.ends          = [];
NB.CHANNELSFIN   = [];
NB.DOMCHANNFIN   = [];
NB.ISINThreshold = nan;
end
