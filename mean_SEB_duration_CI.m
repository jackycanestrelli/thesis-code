function mean_SEB_duration_CI()
% mean_SEB_duration_CI  Mean single-electrode burst duration across DIVs with IQR outlier removal.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, SampFreq, ParamsBurst
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','SampFreq','ParamsBurst'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% MEAN SINGLE-ELECTRODE BURST DURATION ACROSS SELECTED DIVs
% + IQR OUTLIER REMOVAL
% Per ogni well:
%   - split spike per elettrodo
%   - run DetectBurstsNeuro su ogni elettrodo
%   - durata burst = (endSample - startSample) / SampFreq  [s]
%   - meanBurstDurWell = media di tutte le durate di tutti gli elettrodi

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
usePercentileFilter = false;
percentileThr       = 15;
outlierK            = 1.5;      % IQR factor: 1.5 = standard, 3 = solo estremi
recordingTime       = 600;
saveFigure          = true;
saveFolder          = fullfile(outRoot, 'all_outputs', 'SEB_duration');
figName             = 'MeanSingleElectrodeBurstDuration_4CellLines';

if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

if ~exist('ParamsBurst','var')
    error(['ParamsBurst non trovato nel workspace.' newline ...
           'Crea ParamsBurst prima di eseguire questo script.']);
end
if ~exist('SampFreq','var'), error('SampFreq non trovato nel workspace.'); end
if ~exist('recDuration','var')
    warning('recDuration non trovato. Uso recordingTime = %d s.', recordingTime);
    recDuration = recordingTime;
end

%% ===== BUILD CELL LINE LIST =====
CellLines = {};
if exist('DD1C',  'var'), CellLines{end+1} = struct('name','DD1C',  'Data',DD1C);   end
if exist('DD5A',  'var'), CellLines{end+1} = struct('name','DD5A',  'Data',DD5A);   end
if exist('TUBA',  'var'), CellLines{end+1} = struct('name','TUBA',  'Data',TUBA);   end
if exist('C10902','var'), CellLines{end+1} = struct('name','C10902','Data',C10902); end
if isempty(CellLines), error('Nessuna cell line valida trovata nel workspace.'); end

%% ===== PREPARE OUTPUT =====
nCL  = numel(CellLines);
nDIV = numel(selectedDIVs);

MeanDur_single    = nan(nCL, nDIV);
SEM_Dur    = nan(nCL, nDIV);
N_before   = zeros(nCL, nDIV);
N_after    = zeros(nCL, nDIV);
Thresholds = nan(nCL, nDIV);
ExcludedWells = struct();

%% ===== MAIN LOOP =====
for c = 1:nCL
    cellLineName = CellLines{c}.name;
    cellData     = CellLines{c}.Data;
    plateNames   = fieldnames(cellData);

    for d = 1:nDIV
        DIV = selectedDIVs{d};

        dur_values_all = [];
        mfr_values_all = [];
        meta_plate     = {};
        meta_well      = {};

        for p = 1:numel(plateNames)
            meaPlate = plateNames{p};
            if ~isfield(cellData.(meaPlate), DIV), continue; end

            wellNames = fieldnames(cellData.(meaPlate).(DIV));

            for w = 1:numel(wellNames)
                wellName   = wellNames{w};
                spikeTable = cellData.(meaPlate).(DIV).(wellName);

                if ~istable(spikeTable) || height(spikeTable)==0 || ...
                   ~ismember('Channel',spikeTable.Properties.VariableNames) || ...
                   ~ismember('Time',   spikeTable.Properties.VariableNames), continue; end

                Time    = spikeTable.Time;
                Channel = spikeTable.Channel;
                valid   = isfinite(Time) & Time >= 0;
                Time    = Time(valid);
                Channel = Channel(valid);
                if isempty(Time), continue; end

                mfr_well           = numel(Time) / recordingTime;
                wellBurstDurations = getMeanBurstDurationWell(...
                    Time, Channel, SampFreq, ParamsBurst, recDuration);
                if isnan(wellBurstDurations), continue; end

                dur_values_all(end+1) = wellBurstDurations; %#ok<AGROW>
                mfr_values_all(end+1) = mfr_well;           %#ok<AGROW>
                meta_plate{end+1}     = meaPlate;            %#ok<AGROW>
                meta_well{end+1}      = wellName;            %#ok<AGROW>
            end
        end

        N_before(c,d) = numel(dur_values_all);
        if isempty(dur_values_all), continue; end

        % --- filtro percentile su mFR ---
        if usePercentileFilter
            thr      = prctile(mfr_values_all, percentileThr);
            keepMask = mfr_values_all >= thr;
            Thresholds(c,d) = thr;
        else
            thr      = NaN;
            keepMask = true(size(dur_values_all));
        end

        exclIdx  = find(~keepMask);
        exclList = {};
        for k = 1:numel(exclIdx)
            ii = exclIdx(k);
            exclList{end+1} = sprintf('%s | %s', meta_plate{ii}, meta_well{ii}); %#ok<AGROW>
        end
        ExcludedWells.(cellLineName).(DIV) = exclList;

        keptVals = dur_values_all(keepMask);

        % --- rimuovi outliers IQR ---
        keptVals     = removeOutliersIQR(keptVals, outlierK);
        N_after(c,d) = numel(keptVals);

        if ~isempty(keptVals)
            MeanDur_single(c,d) = mean(keptVals, 'omitnan');
            n__ = numel(keptVals);
            if n__ <= 1
                SEM_Dur(c,d) = NaN;
            else
                SEM_Dur(c,d) = tinv(0.975, n__-1) * std(keptVals, 0, 'omitnan') / sqrt(n__);
            end
end
    end
end

%% ===== PLOT =====
colors = [0.70 0.10 0.10;   % DD1C  (rosso)
          0.95 0.55 0.10;   % DD5A  (arancio)
          0.10 0.30 0.70;   % TUBA  (blu)
          0.30 0.60 0.85];  % C10902 (azzurro)
x       = 1:nDIV;

figure('Color','w','Position',[100 100 900 580]);
hold on; ax = gca;

for c = 1:nCL
    y   = MeanDur_single(c,:);
    sem = SEM_Dur(c,:);
    vm  = ~isnan(y) & ~isnan(sem);
    if ~any(~isnan(y)), continue; end

    % --- Shaded CI band (only on points where SEM is valid) ---
    if any(vm)
        xv = x(vm);
        yv = y(vm);
        ev = sem(vm);
        xfill = [xv, fliplr(xv)];
        yfill = [yv+ev, fliplr(yv-ev)];
        fill(ax, xfill, yfill, colors(c,:), 'FaceAlpha', 0.15, ...
             'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % --- Line + filled circle markers, uniform style for all cell lines ---
    valid = ~isnan(y);
    plot(ax, x(valid), y(valid), '-o', ...
        'Color', colors(c,:), 'LineWidth', 1.8, ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(c,:), ...
        'DisplayName', CellLines{c}.name);
end

xlim([0.5 nDIV+0.5]); xticks(x); xticklabels(selectedDIVs);
xlabel('DIV', 'FontSize', 15);
ylabel('Mean single-electrode burst duration (s) ± 95% CI', 'FontSize', 15);
% title(sprintf('Single-electrode burst duration | outliers removed (IQR k=%.1f)', outlierK), ...
%     'FontSize', 12);
legend('Location','best','FontSize',10);
grid on; box on;

%% ===== SUMMARY =====
fprintf('\n===== SEB DURATION SUMMARY (IQR k=%.1f) =====\n', outlierK);
for c = 1:nCL
    cellLineName = CellLines{c}.name;
    fprintf('\n%s\n', cellLineName);
    for d = 1:nDIV
        DIV = selectedDIVs{d};
        fprintf('  %s: kept=%d/%d | mean=%.4f s | 95CI=%.4f\n', ...
            DIV, N_after(c,d), N_before(c,d), MeanDur_single(c,d), SEM_Dur(c,d));
        if isfield(ExcludedWells,cellLineName) && isfield(ExcludedWells.(cellLineName),DIV)
            el = ExcludedWells.(cellLineName).(DIV);
            if ~isempty(el)
                fprintf('  Excluded (pct filter):\n');
                for k = 1:numel(el), fprintf('    - %s\n', el{k}); end
            end
        end
    end
end

if saveFigure
    exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution',300);
    fprintf('\nSalvato: %s\n', fullfile(saveFolder,[figName '.png']));
end

end

function meanDur = getMeanBurstDurationWell(Time, Channel, SampFreq, ParamsBurst, recDuration)
Channel    = string(Channel);
electrodes = unique(Channel);
allDurations = [];
t0 = min(Time);
for e = 1:numel(electrodes)
    elec  = electrodes(e);
    spT   = sort(Time(Channel == elec) - t0);
    sps   = unique(round(spT * SampFreq) + 1);
    sps   = sps(sps >= 1);
    if numel(sps) < 3, continue; end
    bt    = DetectBurstsNeuro(sps, SampFreq, recDuration, ParamsBurst);
    if isempty(bt) || size(bt,1)==0, continue; end
    durs  = (bt(:,2) - bt(:,1)) / SampFreq;
    allDurations = [allDurations; durs(:)]; %#ok<AGROW>
end
if isempty(allDurations), meanDur = NaN;
else, meanDur = mean(allDurations, 'omitnan'); end
end

function clean = removeOutliersIQR(data, k)
% Tukey IQR outlier removal. k=1.5 standard, k=3 solo estremi.
data = data(~isnan(data));
if numel(data) < 4, clean = data; return; end
q1    = prctile(data, 25);
q3    = prctile(data, 75);
iqr_  = q3 - q1;
clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end
