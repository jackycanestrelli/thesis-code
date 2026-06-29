function T = compute_all_metrics_csv(cellData, nbTampere, nbDoorn, SampFreq, ParamsBurst, cfg)
% compute_all_metrics_csv  Compute per-well activity metrics for every
%   cell line / plate / DIV / well and detector, and return them as one table
%   (optionally written to CSV/MAT).
%
% Metrics per row: MFR (mean firing rate, Hz), SBR (single-electrode burst
%   rate, bursts/min), SBD (single-electrode burst duration, s), NBR (network
%   burst rate, NB/min), FragPerNB (mean fragments among fragmented NBs only),
%   FragPrev (fraction of fragmented NBs, [0,1]), NBD (mean NB duration, s),
%   CV_IBI (CV of inter-burst interval).
%
% INPUT
%   cellData    : struct with one field per cell line (e.g. cellData.DD1C),
%                 each a nested struct .(plate).(DIV).(well) = spike table
%   nbTampere   : struct nbTampere.(cellLine) = Tampere NB results (or [])
%   nbDoorn     : struct nbDoorn.(cellLine)   = Doorn NB results (or [])
%   SampFreq    : sampling frequency [Hz]
%   ParamsBurst : parameter struct for DetectBurstsNeuro
%   cfg         : (optional) struct with fields
%                   .celllines       (default {'DD1C','DD5A','TUBA','C10902'})
%                   .selectedDIVs    (default [35 39 41 ... 90])
%                   .recordingTime_s (default 600)
%                   .saveCSV         (default false)
%                   .saveFolder, .csvName  (used when saveCSV is true)
%
% OUTPUT
%   T : table with columns CellLine, Plate, DIV, Well, Detector, MFR, SBR,
%       SBD, NBR, FragPerNB, FragPrev, NBD, CV_IBI

if nargin < 6, cfg = struct(); end
if ~isfield(cfg,'celllines'),       cfg.celllines = {'DD1C','DD5A','TUBA','C10902'}; end
if ~isfield(cfg,'selectedDIVs'),    cfg.selectedDIVs = [35,39,41,46,48,53,55,60,62,67,69,74,76,81,83,88,90]; end
if ~isfield(cfg,'recordingTime_s'), cfg.recordingTime_s = 600; end
if ~isfield(cfg,'saveCSV'),         cfg.saveCSV = false; end
if ~isfield(cfg,'csvName'),         cfg.csvName = 'metrics_per_well.csv'; end

celllines       = cfg.celllines;
selectedDIVs    = cfg.selectedDIVs;
detectors       = {'TampereV1','Doorn'};
recordingTime_s = cfg.recordingTime_s;
recordingTime_m = recordingTime_s / 60;

if cfg.saveCSV && ~exist(cfg.saveFolder,'dir'), mkdir(cfg.saveFolder); end

T = table('Size',[0 13], ...
    'VariableTypes',{'string','string','double','string','string', ...
                     'double','double','double','double','double','double','double','double'}, ...
    'VariableNames',{'CellLine','Plate','DIV','Well','Detector', ...
                     'MFR','SBR','SBD','NBR','FragPerNB','FragPrev','NBD','CV_IBI'});

for ic = 1:numel(celllines)
    cl = celllines{ic};
    if ~isfield(cellData, cl)
        warning('[%s] spike data not found — skipping.', cl);
        continue
    end
    spkData = cellData.(cl);

    nbT = []; nbD = [];
    if ~isempty(nbTampere) && isfield(nbTampere, cl), nbT = nbTampere.(cl); end
    if ~isempty(nbDoorn)   && isfield(nbDoorn,   cl), nbD = nbDoorn.(cl);   end

    plateNames = fieldnames(spkData);

    for ip = 1:numel(plateNames)
        plate = plateNames{ip};

        for id = 1:numel(selectedDIVs)
            divNum   = selectedDIVs(id);
            divField = sprintf('DIV%d', divNum);
            if ~isfield(spkData.(plate), divField)
                alt = sprintf('DIV_%d', divNum);
                if isfield(spkData.(plate), alt), divField = alt;
                else, continue; end
            end

            wellNames = fieldnames(spkData.(plate).(divField));

            for iw = 1:numel(wellNames)
                well       = wellNames{iw};
                spikeTable = spkData.(plate).(divField).(well);

                % MFR — total spikes in the well / recording duration (thesis 3.3.2)
                mfr = NaN;
                if istable(spikeTable) && height(spikeTable) > 0
                    mfr = height(spikeTable) / recordingTime_s;
                end

                % SBR & SBD
                sbr = NaN;  sbd = NaN;
                if istable(spikeTable) && height(spikeTable) > 0 && ...
                        ismember('Channel',spikeTable.Properties.VariableNames) && ...
                        ismember('Time',   spikeTable.Properties.VariableNames)
                    [sbr, sbd] = computeSingleBurstMetrics( ...
                        spikeTable, SampFreq, ParamsBurst, recordingTime_s);
                end

                % NB metrics per detector
                for idet = 1:numel(detectors)
                    det = detectors{idet};
                    nbr = NaN;  fragPerNB = NaN;  fragPrev = NaN;
                    nbd = NaN;  cv_ibi    = NaN;

                    switch det
                        case 'TampereV1', nbStruct = nbT;
                        case 'Doorn',     nbStruct = nbD;
                        otherwise,        nbStruct = [];
                    end

                    if ~isempty(nbStruct) && isfield(nbStruct, plate) && ...
                            isfield(nbStruct.(plate), divField) && ...
                            isfield(nbStruct.(plate).(divField), well)

                        nb = nbStruct.(plate).(divField).(well);

                        if ~isempty(nb)
                            % NBR
                            if isstruct(nb) && isfield(nb,'starts') && ~isempty(nb.starts)
                                nbr = numel(nb.starts) / recordingTime_m;
                            elseif isstruct(nb) && isfield(nb,'num_bursts')
                                nbr = nb.num_bursts / recordingTime_m;
                            elseif isnumeric(nb) && size(nb,2)==2
                                nbr = size(nb,1) / recordingTime_m;
                            end

                            % FragPrev = fraction of NBs with >1 fragment;
                            % FragPerNB = mean fragment count among those only.
                            if isstruct(nb) && isfield(nb,'num_fragments') && ~isempty(nb.num_fragments)
                                nf = nb.num_fragments(:);
                                nf = nf(~isnan(nf));
                                if ~isempty(nf)
                                    % FragPrev: fraction of fragmented NBs (>= 2 fragments, thesis 3.2.5)
                                    fragPrev = sum(nf > 1) / numel(nf);
                                    % FragPerNB: mean fragments over NBs with >= 1 fragment (thesis 3.3.2)
                                    hasFrag = nf >= 1;
                                    if any(hasFrag)
                                        fragPerNB = mean(nf(hasFrag), 'omitnan');
                                    end
                                end
                            end

                            % NBD
                            if isstruct(nb) && isfield(nb,'starts') && ...
                                    isfield(nb,'ends') && ~isempty(nb.starts)
                                nbd = mean(nb.ends(:) - nb.starts(:),'omitnan');
                            elseif isnumeric(nb) && size(nb,2)==2
                                nbd = mean(nb(:,2) - nb(:,1),'omitnan');
                            end

                            % CV_IBI
                            starts = [];
                            if isstruct(nb) && isfield(nb,'starts') && ~isempty(nb.starts)
                                starts = nb.starts(:);
                            elseif isnumeric(nb) && size(nb,2)==2
                                starts = nb(:,1);
                            end
                            if numel(starts) >= 2
                                ibi = diff(sort(starts));
                                mu_ = mean(ibi,'omitnan');
                                sg_ = std(ibi,0,'omitnan');
                                if mu_ > 0, cv_ibi = sg_ / mu_; end
                            end
                        end
                    end

                    T(end+1,:) = {string(cl), string(plate), divNum, ...
                                  string(well), string(det), ...
                                  mfr, sbr, sbd, nbr, fragPerNB, fragPrev, ...
                                  nbd, cv_ibi}; %#ok<AGROW>
                end
            end
        end
    end
end

if cfg.saveCSV
    outPath = fullfile(cfg.saveFolder, cfg.csvName);
    writetable(T, outPath);
    save(fullfile(cfg.saveFolder,'metrics_per_well.mat'), 'T');
    fprintf('CSV saved: %s\n', outPath);
end
end


function [sbr, sbd] = computeSingleBurstMetrics(spikeTable, SampFreq, ParamsBurst, recTime_s)
% Single-electrode burst rate (bursts/min) and mean burst duration (s) for one well.
recTime_m = recTime_s / 60;
Time    = spikeTable.Time;
Channel = string(spikeTable.Channel);
valid   = isfinite(Time) & Time >= 0;
Time    = Time(valid);
Channel = Channel(valid);
if isempty(Time), sbr = NaN; sbd = NaN; return; end

t0         = min(Time);
electrodes = unique(Channel);
allBursts  = 0;
allDurs    = [];

for e = 1:numel(electrodes)
    spT = sort(Time(Channel == electrodes(e)) - t0);
    sps = unique(round(spT * SampFreq) + 1);
    sps = sps(sps >= 1);
    if numel(sps) < 3, continue; end
    bt = DetectBurstsNeuro(sps, SampFreq, recTime_s, ParamsBurst);
    if isempty(bt) || size(bt,1) == 0, continue; end
    allBursts = allBursts + size(bt,1);
    durs = (bt(:,2) - bt(:,1)) / SampFreq;
    allDurs = [allDurs; durs(:)]; %#ok<AGROW>
end

sbr = allBursts / recTime_m;
sbd = mean(allDurs,'omitnan');
if isempty(allDurs), sbd = NaN; end
end
