function NB_out = run_NBdetection_tampere(CellData, Params, detectorFcn)
% run_NBdetection_tampere  Run the Tampere (log ISI-N) network-burst detector
%   over every MEA plate / DIV / well of one cell line and attach a smoothed
%   firing-rate trace to each well result.
%
% INPUT
%   CellData    : nested struct CellData.(meaPlate).(DIV).(well) = spike table
%                 (columns Time [s] and Channel) for one cell line
%   Params      : detector parameter struct (see DetectNetworkBurstsNeuro);
%                 must also contain timeBin, kernelsize, sigma for the FR trace
%   detectorFcn : (optional) handle to the detector, default @DetectNetworkBurstsNeuro.
%                 Pass @DetectNetworkBurstsNeuroV2 for the V2 threshold logic.
%
% OUTPUT
%   NB_out : nested struct NB_out.(meaPlate).(DIV).(well) = NB struct returned
%            by the detector, with extra fields .spikeTimes and .spikeratesmooth

if nargin < 3 || isempty(detectorFcn)
    detectorFcn = @DetectNetworkBurstsNeuro;
end

NB_out    = struct();
meaPlates = fieldnames(CellData);

for p = 1:numel(meaPlates)
    meaPlate = meaPlates{p};

    divNames = fieldnames(CellData.(meaPlate));
    divNames = divNames(startsWith(divNames, 'DIV', 'IgnoreCase', true));

    for d = 1:numel(divNames)
        divName = divNames{d};
        wells   = fieldnames(CellData.(meaPlate).(divName));

        for w = 1:numel(wells)
            wellName   = wells{w};
            spikeTable = CellData.(meaPlate).(divName).(wellName);

            NB = detectorFcn([], spikeTable, Params);

            % Smoothed firing-rate trace (same kernel as the detector), computed
            % for every well regardless of how many NBs were found.
            if ~isempty(spikeTable) && ismember('Time', spikeTable.Properties.VariableNames)
                allSpikeTimes = sort(spikeTable.Time);
                if numel(allSpikeTimes) >= 2
                    dt          = Params.timeBin;
                    tVec        = (allSpikeTimes(1) : dt : allSpikeTimes(end))';
                    spikeCounts = histcounts(allSpikeTimes, [tVec; tVec(end)+dt])';
                    kernel      = gaussian_kernel(Params.kernelsize, Params.sigma);
                    kernel      = kernel / sum(kernel);
                    NB.spikeTimes      = tVec;
                    NB.spikeratesmooth = smooth_signal(spikeCounts / dt, kernel);
                end
            end

            NB_out.(meaPlate).(divName).(wellName) = NB;
        end
    end
end
end
