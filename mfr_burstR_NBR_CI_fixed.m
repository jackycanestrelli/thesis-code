function mfr_burstR_NBR_CI_fixed()
% mfr_burstR_NBR_CI_fixed  Mean firing rate / burst-rate / NB-rate line plots with 95% CI + IQR outlier removal.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Doorn_DD1C, NB_Doorn_DD5A,
%   NB_Doorn_TUBA, NB_Doorn_C10902, SampFreq, ParamsBurst
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Doorn_DD1C','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902','SampFreq','ParamsBurst'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% MEAN FIRING RATE ACROSS MANUALLY CHOSEN DIVs
% WITH PERCENTILE-BASED LOW-ACTIVITY WELL EXCLUSION + IQR OUTLIER REMOVAL
% mFR per well = total spikes / recording time
% Outlier removal: IQR method (Tukey), k settabile

% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV62','DIV67','DIV83','DIV88'};
DIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
recordingTime_s       = 600;      % seconds (recording duration)
usePercentileFilter = false;
percentileThr       = 15;
outlierK            = 1.5;      % IQR factor: 1.5 = standard, 3 = only extremes
saveFigure          = true;
saveFolder          = fullfile(outRoot, 'all_outputs', 'mfr_BR_NBR_CI');
figName             = 'MeanFiringRate_4CellLines';

if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

%% ===== BUILD CELL LINE LIST =====
CellLines = {};
if exist('DD1C',  'var'), CellLines{end+1} = struct('name','DD1C',  'Data',DD1C);   end
if exist('DD5A',  'var'), CellLines{end+1} = struct('name','DD5A',  'Data',DD5A);   end
if exist('TUBA',  'var'), CellLines{end+1} = struct('name','TUBA',  'Data',TUBA);   end
if exist('C10902','var'), CellLines{end+1} = struct('name','C10902','Data',C10902); end
if isempty(CellLines), error('No valid cell line variables found in workspace.'); end

%% ===== PREPARE OUTPUT =====
nCL  = numel(CellLines);
nDIV = numel(selectedDIVs);

MeanMFR    = nan(nCL, nDIV);
SEM_MFR    = nan(nCL, nDIV);
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
                if ~istable(spikeTable) || height(spikeTable)==0, continue; end
                if ~ismember('Time',    spikeTable.Properties.VariableNames), continue; end
                if ~ismember('Channel', spikeTable.Properties.VariableNames), continue; end

                mfr_well = height(spikeTable) / recordingTime_s;
                if isfinite(mfr_well)
                    mfr_values_all(end+1) = mfr_well; %#ok<AGROW>
                    meta_plate{end+1}     = meaPlate;  %#ok<AGROW>
                    meta_well{end+1}      = wellName;  %#ok<AGROW>
                end
            end
        end

        N_before(c,d) = numel(mfr_values_all);
        if isempty(mfr_values_all), continue; end

        if usePercentileFilter
            thr      = prctile(mfr_values_all, percentileThr);
            keepMask = mfr_values_all >= thr;
            Thresholds(c,d) = thr;
        else
            thr      = NaN;
            keepMask = true(size(mfr_values_all));
        end

        exclIdx  = find(~keepMask);
        exclList = {};
        for k = 1:numel(exclIdx)
            ii = exclIdx(k);
            exclList{end+1} = sprintf('%s | %s', meta_plate{ii}, meta_well{ii}); %#ok<AGROW>
        end
        ExcludedWells.(cellLineName).(DIV) = exclList;

        keptVals     = mfr_values_all(keepMask);
        keptVals     = removeOutliersIQR(keptVals, outlierK);
        N_after(c,d) = numel(keptVals);

        if ~isempty(keptVals)
            MeanMFR(c,d) = mean(keptVals);
                            n__ = numel(keptVals);
                if n__ <= 1
                    SEM_MFR(c,d) = NaN;
                else
                    SEM_MFR(c,d) = tinv(0.975, n__-1) * std(keptVals, 0, 'omitnan') / sqrt(n__);
                end
end
    end
end

%% ===== PLOT MFR =====
colors = [0.70 0.10 0.10;   % DD1C  (severe)
          0.95 0.55 0.10;   % DD5A  (mild)
          0.10 0.30 0.70;   % TUBA  (ctrl)
          0.30 0.60 0.85];  % C10902 (ctrl)
markers = {'o','s','d','^'};
x       = 1:nDIV;

figure('Color','w','Position',[100 100 1000 600]);
hold on; ax = gca;
for c = 1:nCL
    y  = MeanMFR(c,:);
    se = SEM_MFR(c,:);
    if ~any(~isnan(y)), continue; end
    vm = ~isnan(y) & ~isnan(se);
    if any(vm)
        xv = x(vm); yv = y(vm); ev = se(vm);
        fill(ax, [xv, fliplr(xv)], [yv+ev, fliplr(max(yv-ev,0))], colors(c,:), ...
             'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
    end
    valid = ~isnan(y);
    plot(ax, x(valid), y(valid), '-o', 'Color',colors(c,:), 'LineWidth',1.8, ...
        'MarkerSize',6, 'MarkerFaceColor',colors(c,:), ...
        'DisplayName', CellLines{c}.name);
end
xlim([0.5 nDIV+0.5]); xticks(x); xticklabels(DIVs);
xlabel('DIV'); ylabel('Mean firing rate (Hz) ± 95% CI');
% title(sprintf('Mean firing rate | IQR outlier removal (k=%.1f)', outlierK));
legend('Location','best'); grid on; box on;

fprintf('\n===== MFR SUMMARY (IQR k=%.1f) =====\n', outlierK);
for c = 1:nCL
    fprintf('\n%s\n', CellLines{c}.name);
    for d = 1:nDIV
        fprintf('  %s: kept=%d/%d | mean=%.4f Hz | 95CI=%.4f\n', ...
            selectedDIVs{d}, N_after(c,d), N_before(c,d), MeanMFR(c,d), SEM_MFR(c,d));
    end
end
if saveFigure
    exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution',300);
end

% % %% ============================================================
% % %% ============================================================
% % % MEAN SINGLE-ELECTRODE BURST RATE ACROSS SELECTED DIVs
% % %% ============================================================
% % 
% % %% ===== USER SETTINGS =====
% selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
% usePercentileFilter = false;
% recordingTime= 10;
% percentileThr       = 15;
% outlierK            = 1.5;
% saveFigure          = true;
% figName             = 'MeanSingleElectrodeBurstCount_4CellLines';
% if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end
% if ~exist('ParamsBurst','var')
%     error(['ParamsBurst not found in workspace.' newline 'Create ParamsBurst before running.']);
% end
% %% ===== BUILD CELL LINE LIST =====
% CellLines = {};
% if exist('DD1C',  'var'), CellLines{end+1} = struct('name','DD1C',  'Data',DD1C);   end
% if exist('DD5A',  'var'), CellLines{end+1} = struct('name','DD5A',  'Data',DD5A);   end
% if exist('TUBA',  'var'), CellLines{end+1} = struct('name','TUBA',  'Data',TUBA);   end
% if exist('C10902','var'), CellLines{end+1} = struct('name','C10902','Data',C10902); end
% if isempty(CellLines), error('No valid cell line variables found in workspace.'); end
% nCL  = numel(CellLines);
% nDIV = numel(selectedDIVs);
% MeanBurstCount = nan(nCL, nDIV);
% SEMBurstCount  = nan(nCL, nDIV);
% N_before       = zeros(nCL, nDIV);
% N_after        = zeros(nCL, nDIV);
% ExcludedWells  = struct();
% for c = 1:nCL
%     cellLineName = CellLines{c}.name;
%     cellData     = CellLines{c}.Data;
%     plateNames   = fieldnames(cellData);
%     for d = 1:nDIV
%         DIV = selectedDIVs{d};
%         burst_values_all = [];
%         mfr_values_all   = [];
%         meta_plate       = {};
%         meta_well        = {};
%         for p = 1:numel(plateNames)
%             meaPlate = plateNames{p};
%             if ~isfield(cellData.(meaPlate), DIV), continue; end
%             wellNames = fieldnames(cellData.(meaPlate).(DIV));
%             for w = 1:numel(wellNames)
%                 wellName   = wellNames{w};
%                 spikeTable = cellData.(meaPlate).(DIV).(wellName);
%                 if ~istable(spikeTable) || height(spikeTable)==0 || ...
%                    ~ismember('Channel',spikeTable.Properties.VariableNames) || ...
%                    ~ismember('Time',   spikeTable.Properties.VariableNames), continue; end
%                 Time      = spikeTable.Time;
%                 validMask = isfinite(Time) & Time >= 0;
%                 Time      = Time(validMask);
%                 if isempty(Time), continue; end
%                 mfr_well          = numel(Time) / recordingTime;
%                 totalBurstsInWell = getTotalBurstsInWell(spikeTable, SampFreq, ParamsBurst, recordingTime);
%                 if isnan(totalBurstsInWell), continue; end
%                 burst_values_all(end+1) = totalBurstsInWell / recordingTime; %#ok<AGROW>
%                 mfr_values_all(end+1)   = mfr_well;                        %#ok<AGROW>
%                 meta_plate{end+1}       = meaPlate;                         %#ok<AGROW>
%                 meta_well{end+1}        = wellName;                         %#ok<AGROW>
%             end
%         end
%         N_before(c,d) = numel(burst_values_all);
%         if isempty(burst_values_all), continue; end
%         if usePercentileFilter
%             thr      = prctile(mfr_values_all, percentileThr);
%             keepMask = mfr_values_all >= thr;
%         else
%             keepMask = true(size(burst_values_all));
%         end
%         exclIdx  = find(~keepMask);
%         exclList = {};
%         for k = 1:numel(exclIdx)
%             ii = exclIdx(k);
%             exclList{end+1} = sprintf('%s | %s', meta_plate{ii}, meta_well{ii}); %#ok<AGROW>
%         end
%         ExcludedWells.(cellLineName).(DIV) = exclList;
%         keptVals     = burst_values_all(keepMask);
%         keptVals     = removeOutliersIQR(keptVals, outlierK);
%         N_after(c,d) = numel(keptVals);
%         if ~isempty(keptVals)
%             MeanBurstCount(c,d) = mean(keptVals);
%                             n__ = numel(keptVals);
%                 if n__ <= 1
%                     SEMBurstCount(c,d) = NaN;
%                 else
%                     SEMBurstCount(c,d) = tinv(0.975, n__-1) * std(keptVals, 0, 'omitnan') / sqrt(n__);
%                 end
% end
%     end
% end
% %% ===== PLOT SEB RATE =====
% figure('Color','w','Position',[100 100 1000 600]);
% hold on; ax = gca;
% for c = 1:nCL
%     y  = MeanBurstCount(c,:);
%     se = SEMBurstCount(c,:);
%     if ~any(~isnan(y)), continue; end
%     vm = ~isnan(y) & ~isnan(se);
%     if any(vm)
%         xv = x(vm); yv = y(vm); ev = se(vm);
%         fill(ax, [xv, fliplr(xv)], [yv+ev, fliplr(max(yv-ev,0))], colors(c,:), ...
%              'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
%     end
%     valid = ~isnan(y);
%     plot(ax, x(valid), y(valid), '-o', 'Color',colors(c,:), 'LineWidth',1.8, ...
%         'MarkerSize',6, 'MarkerFaceColor',colors(c,:), ...
%         'DisplayName', CellLines{c}.name);
% end
% xlim([0.5 nDIV+0.5]); xticks(x); xticklabels(selectedDIVs);
% xlabel('DIV'); ylabel('Mean burst rate (bursts/min) ± 95% CI');
% title(sprintf('Single-electrode burst rate | IQR outlier removal (k=%.1f)', outlierK));
% legend('Location','best'); grid on; box on;
% fprintf('\n===== SEB RATE SUMMARY (IQR k=%.1f) =====\n', outlierK);
% for c = 1:nCL
%     fprintf('\n%s\n', CellLines{c}.name);
%     for d = 1:nDIV
%         fprintf('  %s: kept=%d/%d | mean=%.4f bursts/min | 95CI=%.4f\n', ...
%             selectedDIVs{d}, N_after(c,d), N_before(c,d), ...
%             MeanBurstCount(c,d), SEMBurstCount(c,d));
%     end
% end
% if saveFigure
%     exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution',300);
% end
% %% ============================================================
% %% ============================================================
% % MEAN NETWORK BURST RATE — TAMPERE
% % Legge da NB_Tampere_<CellLine> (già calcolata nel workspace)
% % Struttura: NB_Tampere_DD1C.(meaPlate).(DIV).(wellName).starts
% %% ============================================================
% %% ===== USER SETTINGS =====
% selectedDIVs    = {'DIV53','DIV62','DIV67','DIV83','DIV88'};
% DIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
% recordingTime       = 10;
% usePercentileFilter = false;
% percentileThr       = 15;
% outlierK            = 1.5;
% saveFigure          = true;
% DIVs={'DIV53','DIV60','DIV67','DIV83','DIV90'};
% figName             = 'MeanNetworkBurstCount_4CellLines';
% if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end
% %% ===== BUILD CELL LINE LIST — usa strutture NB pre-calcolate =====
% CellLines_NB_T = {};
% if exist('NB_Tampere_DD1C',  'var'), CellLines_NB_T{end+1} = struct('name','DD1C',  'NB', NB_Tampere_DD1C);   end
% if exist('NB_Tampere_DD5A',  'var'), CellLines_NB_T{end+1} = struct('name','DD5A',  'NB', NB_Tampere_DD5A);   end
% if exist('NB_Tampere_TUBA',  'var'), CellLines_NB_T{end+1} = struct('name','TUBA',  'NB', NB_Tampere_TUBA);   end
% if exist('NB_Tampere_C10902','var'), CellLines_NB_T{end+1} = struct('name','C10902','NB', NB_Tampere_C10902); end
% if isempty(CellLines_NB_T)
%     error('Nessuna struttura NB_Tampere_* trovata nel workspace. Esegui prima il detector Tampere.');
% end
% nCL  = numel(CellLines_NB_T);
% nDIV = numel(selectedDIVs);
% MeanNBCount = nan(nCL, nDIV);
% SEM_NBCount = nan(nCL, nDIV);
% N_before    = zeros(nCL, nDIV);
% N_after     = zeros(nCL, nDIV);
% for c = 1:nCL
%     cellLineName = CellLines_NB_T{c}.name;
%     NB_Cell      = CellLines_NB_T{c}.NB;
%     plateNames   = fieldnames(NB_Cell);
%     for d = 1:nDIV
%         DIV = selectedDIVs{d};
%         nb_values_all = [];
%         for p = 1:numel(plateNames)
%             meaPlate = plateNames{p};
%             if ~isfield(NB_Cell.(meaPlate), DIV), continue; end
%             wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));
%             for w = 1:numel(wellNames)
%                 wellName = wellNames{w};
%                 NB       = NB_Cell.(meaPlate).(DIV).(wellName);
%                 % conta gli NB: usa .starts se disponibile, altrimenti .num_bursts
%                 if isempty(NB), continue; end
%                 if isfield(NB,'starts') && ~isempty(NB.starts)
%                     totalNBs = numel(NB.starts);
%                 elseif isfield(NB,'num_bursts')
%                     totalNBs = NB.num_bursts;
%                 else
%                     continue
%                 end
%                 if totalNBs == 0
%     fprintf('[INACTIVE] %s | %s — 0 NB, excluded.\n', meaPlate, wellName);
%     continue
% end
% nb_values_all(end+1) = totalNBs / recordingTime; %#ok<AGROW>
%             end
%         end
%         N_before(c,d) = numel(nb_values_all);
%         if isempty(nb_values_all), continue; end
%         keptVals     = removeOutliersIQR(nb_values_all, outlierK);
%         N_after(c,d) = numel(keptVals);
%         if ~isempty(keptVals)
%             MeanNBCount(c,d) = mean(keptVals);
%                             n__ = numel(keptVals);
%                 if n__ <= 1
%                     SEM_NBCount(c,d) = NaN;
%                 else
%                     SEM_NBCount(c,d) = tinv(0.975, n__-1) * std(keptVals, 0, 'omitnan') / sqrt(n__);
%                 end
% end
%     end
% end
% %% ===== PLOT NB TAMPERE =====
% colors = [0.70 0.10 0.10;   % DD1C  (severe)
%           0.95 0.55 0.10;   % DD5A  (mild)
%           0.10 0.30 0.70;   % TUBA  (ctrl)
%           0.30 0.60 0.85];  % C10902 (ctrl)
% markers = {'o','s','d','^'};
% x       = 1:nDIV;
% figure('Color','w','Position',[100 100 1000 600]);
% hold on; ax = gca;
% for c = 1:nCL
%     y  = MeanNBCount(c,:);
%     se = SEM_NBCount(c,:);
%     if ~any(~isnan(y)), continue; end
%     vm = ~isnan(y) & ~isnan(se);
%     if any(vm)
%         xv = x(vm); yv = y(vm); ev = se(vm);
%         fill(ax, [xv, fliplr(xv)], [yv+ev, fliplr(max(yv-ev,0))], colors(c,:), ...
%              'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
%     end
%     valid = ~isnan(y);
%     plot(ax, x(valid), y(valid), '-o', 'Color',colors(c,:), 'LineWidth',1.8, ...
%         'MarkerSize',6, 'MarkerFaceColor',colors(c,:), ...
%         'DisplayName', CellLines_NB_T{c}.name);
% end
% xlim([0.5 nDIV+0.5]); xticks(x); xticklabels(DIVs);
% xlabel('DIV'); ylabel('Mean network burst rate (NB/min) ± 95% CI');
% title(sprintf('Network burst rate — Tampere | IQR outlier removal (k=%.1f)', outlierK));
% legend('Location','best'); grid on; box on;
% fprintf('\n===== NB RATE TAMPERE SUMMARY (IQR k=%.1f) =====\n', outlierK);
% for c = 1:nCL
%     fprintf('\n%s\n', CellLines_NB_T{c}.name);
%     for d = 1:nDIV
%         fprintf('  %s: kept=%d/%d | mean=%.4f NB/min | 95CI=%.4f\n', ...
%             selectedDIVs{d}, N_after(c,d), N_before(c,d), ...
%             MeanNBCount(c,d), SEM_NBCount(c,d));
%     end
% end
% if saveFigure
%     exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution',300);
% end
% %% ============================================================
% %% ============================================================
% % MEAN NETWORK BURST RATE — DOORN
% % Legge da NB_Doorn_<CellLine> (già calcolata nel workspace)
% % Struttura: NB_Doorn_DD1C.(meaPlate).(DIV).(wellName)
% %            -> matrice NBs (Nx2) con [start, stop] in secondi
% %            oppure campo .num_bursts
% %% ============================================================
% %% ===== USER SETTINGS =====
% selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
% recordingTime       = 10;
% usePercentileFilter = false;
% percentileThr       = 15;
% outlierK            = 1.5;
% saveFigure          = true;
% figName             = 'Doorn_MeanNetworkBurstCount_4CellLines';
% if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end
% %% ===== BUILD CELL LINE LIST — usa strutture NB pre-calcolate =====
% CellLines_NB_D = {};
% if exist('NB_Doorn_DD1C',  'var'), CellLines_NB_D{end+1} = struct('name','DD1C',  'NB', NB_Doorn_DD1C);   end
% if exist('NB_Doorn_DD5A',  'var'), CellLines_NB_D{end+1} = struct('name','DD5A',  'NB', NB_Doorn_DD5A);   end
% if exist('NB_Doorn_TUBA',  'var'), CellLines_NB_D{end+1} = struct('name','TUBA',  'NB', NB_Doorn_TUBA);   end
% if exist('NB_Doorn_C10902','var'), CellLines_NB_D{end+1} = struct('name','C10902','NB', NB_Doorn_C10902); end
% if isempty(CellLines_NB_D)
%     error('Nessuna struttura NB_Doorn_* trovata nel workspace. Esegui prima il detector Doorn.');
% end
% nCL  = numel(CellLines_NB_D);
% nDIV = numel(selectedDIVs);
% MeanNBCount_Doorn = nan(nCL, nDIV);
% SEM_NBCount_Doorn = nan(nCL, nDIV);
% N_before          = zeros(nCL, nDIV);
% N_after           = zeros(nCL, nDIV);
% for c = 1:nCL
%     cellLineName = CellLines_NB_D{c}.name;
%     NB_Cell      = CellLines_NB_D{c}.NB;
%     plateNames   = fieldnames(NB_Cell);
%     for d = 1:nDIV
%         DIV = selectedDIVs{d};
%         nb_values_all = [];
%         for p = 1:numel(plateNames)
%             meaPlate = plateNames{p};
%             if ~isfield(NB_Cell.(meaPlate), DIV), continue; end
%             wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));
%             for w = 1:numel(wellNames)
%                 wellName = wellNames{w};
%                 NB       = NB_Cell.(meaPlate).(DIV).(wellName);
%                 if isempty(NB), continue; end
%                 % Doorn: struttura con campo .starts (vettore) o matrice Nx2
%                 if isstruct(NB)
%                     if isfield(NB,'starts') && ~isempty(NB.starts)
%                         totalNBs = numel(NB.starts);
%                     elseif isfield(NB,'num_bursts')
%                         totalNBs = NB.num_bursts;
%                     else
%                         continue
%                     end
%                 elseif isnumeric(NB) && size(NB,2) == 2
%                     % matrice Nx2 [start stop]
%                     totalNBs = size(NB,1);
%                 else
%                     continue
%                 end
%                 if totalNBs == 0
%     fprintf('[INACTIVE] %s | %s — 0 NB, excluded.\n', meaPlate, wellName);
%     continue
% end
% nb_values_all(end+1) = totalNBs / recordingTime; %#ok<AGROW>
%             end
%         end
%         N_before(c,d) = numel(nb_values_all);
%         if isempty(nb_values_all), continue; end
%         keptVals     = removeOutliersIQR(nb_values_all, outlierK);
%         N_after(c,d) = numel(keptVals);
%         if ~isempty(keptVals)
%             MeanNBCount_Doorn(c,d) = mean(keptVals);
%                             n__ = numel(keptVals);
%                 if n__ <= 1
%                     SEM_NBCount_Doorn(c,d) = NaN;
%                 else
%                     SEM_NBCount_Doorn(c,d) = tinv(0.975, n__-1) * std(keptVals, 0, 'omitnan') / sqrt(n__);
%                 end
% end
%     end
% end
% %% ===== PLOT NB DOORN =====
% colors = [0.70 0.10 0.10;   % DD1C  (severe)
%           0.95 0.55 0.10;   % DD5A  (mild)
%           0.10 0.30 0.70;   % TUBA  (ctrl)
%           0.30 0.60 0.85];  % C10902 (ctrl)
% markers = {'o','s','d','^'};
% x       = 1:nDIV;
% figure('Color','w','Position',[100 100 1000 600]);
% hold on
% for c = 1:nCL
%     y  = MeanNBCount_Doorn(c,:);
%     se = SEM_NBCount_Doorn(c,:);
%     if ~any(~isnan(y)), continue; end
%     vm = ~isnan(y) & ~isnan(se);
%     if any(vm)
%         xv = x(vm); yv = y(vm); ev = se(vm);
%         fill(gca, [xv, fliplr(xv)], [yv+ev, fliplr(max(yv-ev,0))], colors(c,:), ...
%              'FaceAlpha',0.15, 'EdgeColor','none', 'HandleVisibility','off');
%     end
%     valid = ~isnan(y);
%     plot(gca, x(valid), y(valid), '-o', 'Color',colors(c,:), 'LineWidth',1.8, ...
%         'MarkerSize',6, 'MarkerFaceColor',colors(c,:), ...
%         'DisplayName', CellLines_NB_D{c}.name);
% end
% xlim([0.5 nDIV+0.5]); xticks(x); xticklabels(DIVs);
% xlabel('DIV'); ylabel('Mean network burst rate (NB/min) ± 95% CI');
% title(sprintf('Network burst rate — Doorn | IQR outlier removal (k=%.1f)', outlierK));
% legend('Location','best'); grid on; box on;
% fprintf('\n===== NB RATE DOORN SUMMARY (IQR k=%.1f) =====\n', outlierK);
% for c = 1:nCL
%     fprintf('\n%s\n', CellLines_NB_D{c}.name);
%     for d = 1:nDIV
%         fprintf('  %s: kept=%d/%d | mean=%.4f NB/min | 95CI=%.4f\n', ...
%             selectedDIVs{d}, N_after(c,d), N_before(c,d), ...
%             MeanNBCount_Doorn(c,d), SEM_NBCount_Doorn(c,d));
%     end
% end
% if saveFigure
%     exportgraphics(gcf, fullfile(saveFolder,[figName '.png']), 'Resolution',300);
% end

% LOCAL FUNCTIONS

end

function totalBurstsInWell = getTotalBurstsInWell(spikeTable, SampFreq, ParamsBurst, recordingTime)
if ~istable(spikeTable) || height(spikeTable)==0 || ...
   ~ismember('Channel',spikeTable.Properties.VariableNames) || ...
   ~ismember('Time',   spikeTable.Properties.VariableNames)
    totalBurstsInWell = NaN; return
end
Time    = spikeTable.Time;
Channel = string(spikeTable.Channel);
valid   = isfinite(Time) & Time >= 0;
Time    = Time(valid);  Channel = Channel(valid);
if isempty(Time), totalBurstsInWell = 0; return; end
t0         = min(Time);
electrodes = unique(Channel);
totalBurstsInWell = 0;
for e = 1:numel(electrodes)
    spT = sort(Time(Channel==electrodes(e)) - t0);
    sps = unique(round(spT * SampFreq) + 1);
    sps = sps(sps >= 1);
    if numel(sps) < 3, continue; end
    bt  = DetectBurstsNeuro(sps, SampFreq, recordingTime, ParamsBurst);
    totalBurstsInWell = totalBurstsInWell + size(bt,1);
end
end

function clean = removeOutliersIQR(data, k)
data  = data(~isnan(data));
if numel(data) < 4, clean = data; return; end
q1    = prctile(data, 25);
q3    = prctile(data, 75);
iqr_  = q3 - q1;
clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end
