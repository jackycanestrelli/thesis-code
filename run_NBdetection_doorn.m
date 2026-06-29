function NB_out = run_NBdetection_doorn(CellData, P)
% run_NBdetection_doorn  Run the Doorn (2024) network-burst detector over every
%   MEA plate / DIV / well of one cell line: bins spikes, builds firing-rate
%   traces, detects network bursts and post-processes them into an NB struct.
%
% INPUT
%   CellData : nested struct CellData.(meaPlate).(DIV).(well) = spike table
%              (columns Time [s] and Channel) for one cell line
%   P        : parameter struct with fields
%                .maxTime          recording stop time [s]
%                .timeBin          bin width [s]
%                .movingWindow_sec moving window for local max [s]
%                .sigma            Gaussian smoothing std [samples]
%                .kernelsize       Gaussian kernel width [samples]
%
% OUTPUT
%   NB_out : nested struct NB_out.(meaPlate).(DIV).(well) = NB_D struct with
%            fields .starts/.ends/.dur [s], .spikes, .num_fragments,
%            .is_fragmented, .fragmentPeakTimes, peak/rate traces (*_Doorn)

NB_out    = struct();
meaPlates = fieldnames(CellData);

for p = 1:numel(meaPlates)
    meaPlate = meaPlates{p};

    divNames = fieldnames(CellData.(meaPlate));
    divNames = divNames(startsWith(divNames, 'DIV', 'IgnoreCase', true));

    for d = 1:numel(divNames)
        DIV   = divNames{d};
        wells = fieldnames(CellData.(meaPlate).(DIV));

        for w = 1:numel(wells)
            wellName   = wells{w};
            spikesData = CellData.(meaPlate).(DIV).(wellName);

            % 1) Binning
            Channel = spikesData.Channel;
            Time    = spikesData.Time;

            electrodes  = unique(Channel);
            nElectrodes = numel(electrodes);

            startTime = 0;
            stopTime  = P.maxTime;
            rectime   = stopTime - startTime;

            validMask   = (Time > startTime) & (Time < stopTime);
            Time_wot    = Time(validMask) - startTime;
            Channel_wot = Channel(validMask);

            nBins     = floor(rectime / P.timeBin);
            APsinbin  = zeros(nElectrodes, nBins);
            spikerate = zeros(nElectrodes, nBins);

            for e = 1:nElectrodes
                APt = Time_wot(strcmp(Channel_wot, electrodes{e}));
                for k = 1:nBins
                    t0 = (k-1)*P.timeBin;
                    t1 =  k   *P.timeBin;
                    APsinbin(e,k)  = sum((APt > t0) & (APt < t1));
                    spikerate(e,k) = APsinbin(e,k) / P.timeBin;
                end
            end

            % 2) Firing rate (weighted traces are used only for plotting)
            APsinbintot          = sum(APsinbin,  1);
            spikeratetot         = sum(spikerate, 1);
            activeChannelsPerBin = sum(APsinbin > 0, 1);
            spikeratetot_weighted = spikeratetot .* activeChannelsPerBin;

            windowBins_plot        = ceil(0.1 / P.timeBin);
            spikerateplot          = movmean(spikeratetot, windowBins_plot);
            spikerateplot_weighted = movmean(spikeratetot_weighted, windowBins_plot);

            kernel = gaussian_kernel(P.kernelsize, P.sigma);
            kernel = kernel / sum(kernel);
            spikeratesmooth_base     = smooth_signal(spikeratetot, kernel);
            spikeratesmooth_weighted = smooth_signal(spikeratetot_weighted, kernel);

            % 3) Doorn detection
            [NBs_doorn, peaksdoorn, spikeratesmooth_detect] = ...
                detectNetworkBursts_Doorn2024(spikeratetot, APsinbin, ...
                P.timeBin, P.movingWindow_sec, P.sigma, P.kernelsize);

            % 4) Post-process (convert bin indices to times; 0-based axis)
            NBs0 = NBs_doorn;
            NBs0(:,1:2)  = NBs0(:,1:2) - 1;
            peakTimesAll = (peaksdoorn - 1) * P.timeBin;

            NB_D = struct();
            NB_D.starts        = NBs0(:,1) * P.timeBin;
            NB_D.ends          = NBs0(:,2) * P.timeBin;
            NB_D.dur           = (NBs0(:,2)-NBs0(:,1)) * P.timeBin;
            NB_D.spikes        = NBs_doorn(:,3);
            NB_D.num_fragments = NBs_doorn(:,5);

            NB_D.spikeratetot_frag_Doorn      = spikeratetot;
            NB_D.spikeratesmooth_Doorn        = spikeratesmooth_base;
            NB_D.spikerateplot_Doorn          = spikerateplot;
            NB_D.spikeratesmooth_detect_Doorn = spikeratesmooth_detect;
            NB_D.spikeTimes_Doorn             = (0:length(spikeratetot)-1) * P.timeBin;
            NB_D.peakTimes_Doorn              = peakTimesAll;

            NB_D.activeChannelsPerBin_Doorn     = activeChannelsPerBin;
            NB_D.spikeratetot_weighted_Doorn    = spikeratetot_weighted;
            NB_D.spikerateplot_weighted_Doorn   = spikerateplot_weighted;
            NB_D.spikeratesmooth_weighted_Doorn = spikeratesmooth_weighted;
            NB_D.APsinbintot_Doorn              = APsinbintot;

            NB_D.fragmentPeakTimes = cell(size(NBs0,1),1);
            for iNB = 1:size(NBs0,1)
                mask = (peakTimesAll >= NB_D.starts(iNB)) & (peakTimesAll < NB_D.ends(iNB));
                NB_D.fragmentPeakTimes{iNB} = peakTimesAll(mask);
            end
            NB_D.is_fragmented = NB_D.num_fragments >= 2;

            NB_out.(meaPlate).(DIV).(wellName) = NB_D;
        end
    end
end
end
