function mean_perc_rnd_spikes_fixed()
% mean_perc_rnd_spikes_fixed  Mean percentage of random (out-of-burst) spikes across DIVs (boxplot).
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

% MEAN PERCENTAGE OF RANDOM SPIKES ACROSS SELECTED DIVs
% + IQR OUTLIER REMOVAL
% + BOXPLOT
% Per ogni well:
% - split spike per elettrodo
% - run DetectBurstsNeuro su ogni elettrodo
% - conta spike totali e spike che cadono dentro burst
% - % random spikes = ((tot spikes - spikes in burst) / tot spikes) * 100
% Output:
% - boxplot dei valori per well per ogni cell line e DIV

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV88'};
DIVs = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
usePercentileFilter = false;
percentileThr = 15;
outlierK = 1.5;  % IQR factor
recordingTime = 600;
saveFigure = true;
saveFolder = fullfile(outRoot, 'all_outputs', 'perc_rnd_spikes');
figName = 'MeanPercentageRandomSpikes_Boxplot_4CellLines';

if saveFigure && ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

if ~exist('ParamsBurst','var')
    error(['ParamsBurst non trovato nel workspace.' newline ...
        'Crea ParamsBurst prima di eseguire questo script.']);
end
if ~exist('SampFreq','var')
    error('SampFreq non trovato nel workspace.');
end
if ~exist('recDuration','var')
    warning('recDuration non trovato. Uso recordingTime = %d s.', recordingTime);
    recDuration = recordingTime;
end

%% ===== BUILD CELL LINE LIST =====
CellLines = {};
if exist('DD1C', 'var'),   CellLines{end+1} = struct('name','DD1C',   'Data',DD1C); end
if exist('DD5A', 'var'),   CellLines{end+1} = struct('name','DD5A',   'Data',DD5A); end
if exist('TUBA', 'var'),   CellLines{end+1} = struct('name','TUBA',   'Data',TUBA); end
if exist('C10902','var'),  CellLines{end+1} = struct('name','C10902', 'Data',C10902); end
if isempty(CellLines)
    error('Nessuna cell line valida trovata nel workspace.');
end

%% ===== PREPARE OUTPUT =====
nCL = numel(CellLines);
nDIV = numel(selectedDIVs);

RndSpikeValues = cell(nCL, nDIV);
N_before = zeros(nCL, nDIV);
N_after = zeros(nCL, nDIV);
Thresholds = nan(nCL, nDIV);
ExcludedWells = struct();

%% ===== MAIN LOOP =====
for c = 1:nCL
    cellLineName = CellLines{c}.name;
    cellData = CellLines{c}.Data;
    plateNames = fieldnames(cellData);

    for d = 1:nDIV
        DIV = selectedDIVs{d};

        rnd_values_all = [];
        mfr_values_all = [];
        meta_plate = {};
        meta_well = {};

        for p = 1:numel(plateNames)
            meaPlate = plateNames{p};
            if ~isfield(cellData.(meaPlate), DIV), continue; end

            wellNames = fieldnames(cellData.(meaPlate).(DIV));

            for w = 1:numel(wellNames)
                wellName = wellNames{w};
                spikeTable = cellData.(meaPlate).(DIV).(wellName);

                if ~istable(spikeTable) || height(spikeTable)==0 || ...
                        ~ismember('Channel', spikeTable.Properties.VariableNames) || ...
                        ~ismember('Time', spikeTable.Properties.VariableNames)
                    continue;
                end

                Time = spikeTable.Time;
                Channel = spikeTable.Channel;
                valid = isfinite(Time) & Time >= 0;
                Time = Time(valid);
                Channel = Channel(valid);

                if isempty(Time), continue; end

                mfr_well = numel(Time) / recordingTime;
                pctRandomWell = getMeanRandomSpikePercentageWell( ...
                    Time, Channel, SampFreq, ParamsBurst, recDuration);

                if isnan(pctRandomWell), continue; end

                rnd_values_all(end+1) = pctRandomWell; %#ok<AGROW>
                mfr_values_all(end+1) = mfr_well; %#ok<AGROW>
                meta_plate{end+1} = meaPlate; %#ok<AGROW>
                meta_well{end+1} = wellName; %#ok<AGROW>
            end
        end

        N_before(c,d) = numel(rnd_values_all);
        if isempty(rnd_values_all), continue; end

        % --- filtro percentile su mFR ---
        if usePercentileFilter
            thr = prctile(mfr_values_all, percentileThr);
            keepMask = mfr_values_all >= thr;
            Thresholds(c,d) = thr;
        else
            thr = NaN;
            keepMask = true(size(rnd_values_all));
        end

        exclIdx = find(~keepMask);
        exclList = {};
        for k = 1:numel(exclIdx)
            ii = exclIdx(k);
            exclList{end+1} = sprintf('%s | %s', meta_plate{ii}, meta_well{ii}); %#ok<AGROW>
        end
        ExcludedWells.(cellLineName).(DIV) = exclList;

        keptVals = rnd_values_all(keepMask);

        % --- rimuovi outliers IQR ---
        keptVals = removeOutliersIQR(keptVals, outlierK);
        N_after(c,d) = numel(keptVals);

        RndSpikeValues{c,d} = keptVals;
    end
end

%% ===== PLOT BOXPLOT =====
colors = [0.70 0.10 0.10;   % DD1C  (rosso)
          0.95 0.55 0.10;   % DD5A  (arancio)
          0.10 0.30 0.70;   % TUBA  (blu)
          0.30 0.60 0.85];  % C10902 (azzurro)

figure('Color','w','Position',[100 100 1400 600]);

for d = 1:nDIV
    subplot(1, nDIV, d); hold on;
% --- fix overlap sgtitle/title: abbassa i subplot lasciando margine in cima ---
ax_pos = get(gca, 'Position');
top_margin = 0.13;   % spazio riservato in cima per sgtitle + title subplot
ax_new_h = ax_pos(4) - top_margin + ax_pos(2) - 0.02;
if ax_new_h > 0.1
    set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3), ax_new_h]);
end

    dataVec = [];
    groupVec = [];
    colorVec = [];

    for c = 1:nCL
        vals = RndSpikeValues{c,d};
        if isempty(vals), continue; end

        nvals = numel(vals);
        dataVec = [dataVec; vals(:)];
        groupVec = [groupVec; repmat(c, nvals, 1)];
        colorVec = [colorVec; repmat(colors(c,:), nvals, 1)];
    end

    if isempty(dataVec)
%         title(DIVs{d});
        ylabel('% random spikes');
        xticks(1:nCL);
        xticklabels(cellfun(@(x)x.name, CellLines, 'UniformOutput', false));
        box on;
        continue;
    end

    boxplot(dataVec, groupVec, ...
        'Labels', cellfun(@(x)x.name, CellLines, 'UniformOutput', false), ...
        'Symbol', 'k+', ...
        'Whisker', 1.5);
       ylim([0 100]);

for c = 1:nCL
    vals = RndSpikeValues{c,d};
    if isempty(vals), continue; end

    m = mean(vals,'omitnan');
    q1 = prctile(vals,25);
    q3 = prctile(vals,75);

    dy = 3;  % percentage points
    if q3 + dy <= 98
        yText = q3 + dy;
        vAlign = 'bottom';
    else
        yText = q1 - dy;
        vAlign = 'top';
    end

    text(c, yText, sprintf('%.2f', m), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment',vAlign, ...
        'FontSize',9, ...
        'FontWeight','bold', ...
        'Color',[0.2 0.2 0.2], ...
        'BackgroundColor','w', ...
        'Margin',1, ...
        'Clipping','on');
end
    h = findobj(gca,'Tag','Box');
    for j = 1:min(numel(h), nCL)
        patch(get(h(j),'XData'), get(h(j),'YData'), colors(numel(h)-j+1,:), ...
            'FaceAlpha', 0.35, 'EdgeColor', colors(numel(h)-j+1,:), 'LineWidth', 1.2);
    end

%     title(DIVs{d}, 'FontSize', 12);
    ylabel('% random spikes', 'FontSize', 15);
    ylim([0 100]);
    grid on;
    box on;
end

% sgtitle(sprintf('Mean percentage of random spikes | outliers removed (IQR k=%.1f)', outlierK), ...
%     'FontSize', 13);

%% ===== SUMMARY =====
fprintf('\n===== RANDOM SPIKES SUMMARY (IQR k=%.1f) =====\n', outlierK);
for c = 1:nCL
    cellLineName = CellLines{c}.name;
    fprintf('\n%s\n', cellLineName);

    for d = 1:nDIV
        DIV = selectedDIVs{d};
        vals = RndSpikeValues{c,d};

        if isempty(vals)
            fprintf('  %s: kept=%d/%d | no valid wells\n', DIV, N_after(c,d), N_before(c,d));
        else
            fprintf('  %s: kept=%d/%d | mean=%.2f %% | median=%.2f %%\n', ...
                DIV, N_after(c,d), N_before(c,d), mean(vals,'omitnan'), median(vals,'omitnan'));
        end

        if isfield(ExcludedWells,cellLineName) && isfield(ExcludedWells.(cellLineName),DIV)
            el = ExcludedWells.(cellLineName).(DIV);
            if ~isempty(el)
                fprintf('    Excluded (pct filter):\n');
                for k = 1:numel(el)
                    fprintf('    - %s\n', el{k});
                end
            end
        end
    end
end

if saveFigure
    exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution', 300);
    fprintf('\nSalvato: %s\n', fullfile(saveFolder,[figName '.png']));
end

end

function pctRandom = getMeanRandomSpikePercentageWell(Time, Channel, SampFreq, ParamsBurst, recDuration)

    Channel = string(Channel);
    electrodes = unique(Channel);

    totalSpikesAll = 0;
    burstSpikesAll = 0;

    t0 = min(Time);

    for e = 1:numel(electrodes)
        elec = electrodes(e);
        spT = sort(Time(Channel == elec) - t0);
        sps = unique(round(spT * SampFreq) + 1);
        sps = sps(sps >= 1);

        if isempty(sps)
            continue;
        end

        totalSpikesAll = totalSpikesAll + numel(sps);

        if numel(sps) < 3
            continue;
        end

        bt = DetectBurstsNeuro(sps, SampFreq, recDuration, ParamsBurst);
        if isempty(bt) || size(bt,1)==0
            continue;
        end

        inBurst = false(size(sps));
        for b = 1:size(bt,1)
            inBurst = inBurst | (sps >= bt(b,1) & sps <= bt(b,2));
        end

        burstSpikesAll = burstSpikesAll + sum(inBurst);
    end

    if totalSpikesAll == 0
        pctRandom = NaN;
    else
        randomSpikes = totalSpikesAll - burstSpikesAll;
        pctRandom = (randomSpikes / totalSpikesAll) * 100;
    end
end

function clean = removeOutliersIQR(data, k)
    data = data(~isnan(data));
    if numel(data) < 4
        clean = data;
        return;
    end
    q1 = prctile(data, 25);
    q3 = prctile(data, 75);
    iqr_ = q3 - q1;
    clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end
