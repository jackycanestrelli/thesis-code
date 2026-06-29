function FR_raster_TvD_modified()
% FR_raster_TvD_modified  Per-well firing-rate + raster figures (Tampere and Doorn side by side).
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Doorn_DD1C, NB_Doorn_DD5A,
%   NB_Doorn_TUBA, NB_Doorn_C10902, maxTime, sigma
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Doorn_DD1C','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902','maxTime','sigma'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% 1 FIGURA PER OGNI WELL (Tampere + Doorn affiancati)
% Per ogni cell line > plate > DIV > well:
%   - 1 figura con 2 pannelli: [Tampere | Doorn]
%   - FR normalizzato con Z-SCORE GLOBALE
%   - Salvata in: Graphs / CellLineName / DIV /
% Struttura dati attesa:
%   DD1C.(meaPlate).(DIV).(well) = table(Channel, Time)
%   NB_Tampere_DD1C.(meaPlate).(DIV).(well)
%   NB_Doorn_DD1C.(meaPlate).(DIV).(well)
%   ... stesso schema per DD5A, TUBA, C10902

%% ===== USER SETTINGS =====
maxTime  = 120;   % secondi — nuovo limite superiore della finestra
tStart   = 60;    % secondi — nuovo limite inferiore della finestra  ← AGGIUNTO
allowedDIVs = [83]; % only DIV83 (all wells, plotted with both Tampere and Doorn)
saveFigures = true;
saveFolder  = fullfile(outRoot, 'all_outputs', 'raster_FR');

set(0,'DefaultFigureVisible','on');
set(0,'DefaultFigureWindowStyle','normal');

%% ===== BUILD CELL LINE CONFIG =====
CL = struct([]);
k  = 0;

if evalin('base','exist(''DD1C'',''var'')') && ...
   evalin('base','exist(''NB_Tampere_DD1C'',''var'')') && ...
   evalin('base','exist(''NB_Doorn_DD1C'',''var'')')
    k = k+1; CL(k).name = 'DD1C';
    CL(k).Data = evalin('base','DD1C');
    CL(k).NB_T = evalin('base','NB_Tampere_DD1C');
    CL(k).NB_D = evalin('base','NB_Doorn_DD1C');
end

if evalin('base','exist(''DD5A'',''var'')') && ...
   evalin('base','exist(''NB_Tampere_DD5A'',''var'')') && ...
   evalin('base','exist(''NB_Doorn_DD5A'',''var'')')
    k = k+1; CL(k).name = 'DD5A';
    CL(k).Data = evalin('base','DD5A');
    CL(k).NB_T = evalin('base','NB_Tampere_DD5A');
    CL(k).NB_D = evalin('base','NB_Doorn_DD5A');
end

if evalin('base','exist(''TUBA'',''var'')') && ...
   evalin('base','exist(''NB_Tampere_TUBA'',''var'')') && ...
   evalin('base','exist(''NB_Doorn_TUBA'',''var'')')
    k = k+1; CL(k).name = 'TUBA';
    CL(k).Data = evalin('base','TUBA');
    CL(k).NB_T = evalin('base','NB_Tampere_TUBA');
    CL(k).NB_D = evalin('base','NB_Doorn_TUBA');
end

if evalin('base','exist(''C10902'',''var'')') && ...
   evalin('base','exist(''NB_Tampere_C10902'',''var'')') && ...
   evalin('base','exist(''NB_Doorn_C10902'',''var'')')
    k = k+1; CL(k).name = 'C10902';
    CL(k).Data = evalin('base','C10902');
    CL(k).NB_T = evalin('base','NB_Tampere_C10902');
    CL(k).NB_D = evalin('base','NB_Doorn_C10902');
end

if isempty(CL)
    error('Nessuna cell line valida trovata nel workspace.');
end

%% ===== FILTRA DIV =====
allDIVs = {};
for c = 1:numel(CL)
    plateNames = fieldnames(CL(c).Data);
    for p = 1:numel(plateNames)
        divNames = fieldnames(CL(c).Data.(plateNames{p}));
        allDIVs  = union(allDIVs, divNames);
    end
end
if isempty(allDIVs), error('Nessun DIV trovato.'); end

allDIVs  = sortDIVs(allDIVs);
keepMask = false(size(allDIVs));
for i = 1:numel(allDIVs)
    keepMask(i) = ismember(getDIVnumber(allDIVs{i}), allowedDIVs);
end
DIVs = allDIVs(keepMask);
if isempty(DIVs), error('Nessun DIV trovato tra quelli specificati in allowedDIVs.'); end

%% ===== CALCOLO Z-SCORE GLOBALE =====
fprintf('Raccolta valori FR per z-score globale...\n');
allFR_vals = [];

for c = 1:numel(CL)
    plateNames = fieldnames(CL(c).Data);
    for p = 1:numel(plateNames)
        meaPlate = plateNames{p};
        if ~isfield(CL(c).Data.(meaPlate), 'dummy')
            divNames = fieldnames(CL(c).Data.(meaPlate));
        else
            continue
        end
        for dv = 1:numel(divNames)
            DIV = divNames{dv};
            if ~ismember(getDIVnumber(DIV), allowedDIVs), continue; end
            wellNames = fieldnames(CL(c).Data.(meaPlate).(DIV));
            for w = 1:numel(wellNames)
                wellName = wellNames{w};
                for mi = 1:2
                    if mi == 1
                        NBroot = CL(c).NB_T; mName = 'Tampere';
                    else
                        NBroot = CL(c).NB_D; mName = 'Doorn';
                    end
                    NB = getNBstruct(NBroot, meaPlate, DIV, wellName);
                    [~, yFR, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
                        extractMethodFields(NB, mName, maxTime);
                    if ~isempty(yFR)
                        allFR_vals = [allFR_vals; yFR(:)]; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

if numel(allFR_vals) < 2
    warning('Pochi valori FR: uso mu=0, sigma=1.');
    global_mu = 0; global_sig = 1;
else
    global_mu  = mean(allFR_vals,'omitnan');
    global_sig = std(allFR_vals, 0,'omitnan');
    if global_sig == 0, global_sig = 1; end
end
fprintf('Z-score globale: mu=%.4f Hz, sigma=%.4f Hz (N=%d)\n', ...
    global_mu, global_sig, numel(allFR_vals));

%% ===== LOOP PRINCIPALE: ogni CL > plate > DIV > well =====
for c = 1:numel(CL)
    cellLineName = CL(c).name;
    plateNames   = fieldnames(CL(c).Data);

    for p = 1:numel(plateNames)
        meaPlate = plateNames{p};
        divNames = fieldnames(CL(c).Data.(meaPlate));

        for dv = 1:numel(divNames)
            DIV = divNames{dv};
            if ~ismember(getDIVnumber(DIV), allowedDIVs), continue; end

            wellNames = fieldnames(CL(c).Data.(meaPlate).(DIV));

            % --- cartella: saveFolder/CellLineName/DIV/
            outDir = fullfile(saveFolder, cellLineName, DIV);
            if saveFigures && ~exist(outDir,'dir')
                mkdir(outDir);
            end

            for w = 1:numel(wellNames)
                wellName   = wellNames{w};
                spikesData = CL(c).Data.(meaPlate).(DIV).(wellName);

                if isempty(spikesData) || height(spikesData) == 0, continue; end
                if ~ismember('Time',spikesData.Properties.VariableNames) || ...
                   ~ismember('Channel',spikesData.Properties.VariableNames), continue; end

                NB_T = getNBstruct(CL(c).NB_T, meaPlate, DIV, wellName);
                NB_D = getNBstruct(CL(c).NB_D, meaPlate, DIV, wellName);

                % --- SALTA well se Tampere O Doorn non hanno NB nella finestra ---
                [~,~,~,~,~, nNB_T_check,~,~,~,~,~] = extractMethodFields(NB_T, 'Tampere', maxTime);
                [~,~,~,~,~, nNB_D_check,~,~,~,~,~] = extractMethodFields(NB_D, 'Doorn',   maxTime);
                if nNB_T_check == 0 || nNB_D_check == 0
                    fprintf('  Skipped (no NB in window): %s/%s/%s/%s  [T=%d, D=%d]\n', ...
                        cellLineName, meaPlate, DIV, wellName, nNB_T_check, nNB_D_check);
                    continue;
                end

                fig = figure( ...
                    'Visible',     'on', ...
                    'WindowStyle', 'normal', ...
                    'Color',       'w', ...
                    'Name',        sprintf('%s | %s | %s | %s', ...
                                   cellLineName, meaPlate, DIV, wellName), ...
                    'NumberTitle', 'off', ...
                    'Position',    [80 80 1400 520]);

                tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

                % --- pannello Tampere ---
                ax1 = nexttile; hold(ax1,'on');
                plotSinglePanel(ax1, spikesData, NB_T, 'Tampere', ...
                    cellLineName, meaPlate, wellName, maxTime,tStart, ...
                    global_mu, global_sig);

                % --- pannello Doorn ---
                ax2 = nexttile; hold(ax2,'on');
                plotSinglePanel(ax2, spikesData, NB_D, 'Doorn', ...
                    cellLineName, meaPlate, wellName, maxTime,tStart, ...
                    global_mu, global_sig);

%                 sgtitle(tl, sprintf('%s | Plate: %s | %s | Well: %s | FR z-scored (global)', ...
%                     cellLineName, meaPlate, DIV, wellName), 'Interpreter','none');

                drawnow;

                if saveFigures
                    fileOut = sprintf('%s_%s_%s_%s_zFR.png', ...
                        cellLineName, meaPlate, DIV, wellName);
                    fileOut = regexprep(fileOut,'[^\w\.-]','_');
                    exportgraphics(fig, fullfile(outDir, fileOut), 'Resolution', 300);
                end

                close(fig);
                fprintf('  Saved: %s/%s/%s/%s\n', cellLineName, meaPlate, DIV, wellName);
            end
        end
    end
end

fprintf('\nDone.\n');

end

function plotSinglePanel(ax, spikesData, NB, methodName, cellLineName, ...
    meaPlate, wellName, maxTime, tStart, global_mu, global_sig)
Channel  = spikesData.Channel;
Time     = spikesData.Time;
maskTime = Time >= tStart & Time <= maxTime;
Channel  = Channel(maskTime);
Time     = Time(maskTime);

if isempty(Time)
    axis(ax,'off');
%     title(ax, sprintf('%s | %s\nno spikes in [0, %ds]', ...
%         methodName, wellName, maxTime), 'Interpreter','none');
    return
end

electrodes  = unique(Channel,'stable');
electrodes  = sortElectrodes(electrodes);
nElectrodes = numel(electrodes);

% --- RASTER sinistro ---
yyaxis(ax,'left')
for e = 1:nElectrodes
    elecSpikes = Time(strcmp(Channel, electrodes{e}));
    plot(ax, elecSpikes, e*ones(numel(elecSpikes),1), ...
        '|', 'MarkerSize', 8, 'Color', [0 0 0]);
end
ylim(ax, [0 nElectrodes+1]);
xlim(ax, [tStart maxTime]);
ylabel(ax, 'Electrodes');

% --- FR z-scored destro ---
yyaxis(ax,'right')
ylabel(ax, 'FR (z-score)');

[tFR_line, yFR_line, tFR_peaks, yFR_peaks, fragPeaks, ...
    nNB, meanFrag, startColor, endColor, ~, peakSpec] = ...
    extractMethodFields(NB, methodName, maxTime);

if ~isempty(yFR_line)
    yFR_line_z  = (yFR_line  - global_mu) / global_sig;
else
    yFR_line_z  = [];
end
if ~isempty(yFR_peaks)
    yFR_peaks_z = (yFR_peaks - global_mu) / global_sig;
else
    yFR_peaks_z = [];
end

if ~isempty(tFR_line) && ~isempty(yFR_line_z)
    plot(ax, tFR_line, yFR_line_z, ...
        'Color', [0.75 0.85 1.00], 'LineWidth', 1.8, 'LineStyle', '--');
    ymax_z = max(yFR_line_z);
    ymin_z = min(yFR_line_z);
    if ymax_z > ymin_z
        ylim(ax, [ymin_z - 0.1*abs(ymin_z+eps), ymax_z*1.1]);
    end
end

if ~isempty(fragPeaks) && ~isempty(tFR_peaks) && ~isempty(yFR_peaks_z)
    peak_vals_z = interp1(tFR_peaks, yFR_peaks_z, fragPeaks, 'linear','extrap');
    if strcmp(peakSpec.type,'filled')
        scatter(ax, fragPeaks, peak_vals_z, peakSpec.size, peakSpec.color, 'filled');
    else
        scatter(ax, fragPeaks, peak_vals_z, peakSpec.size, ...
            'Marker', peakSpec.marker, 'MarkerEdgeColor', peakSpec.color);
    end
end

if ~isempty(NB) && isfield(NB,'starts') && isfield(NB,'ends')
    yyaxis(ax,'left')
    for iNB = 1:numel(NB.starts)
        x1 = NB.starts(iNB);
        x2 = NB.ends(iNB);

        if x2 < 0 || x1 > maxTime
            continue;
        end

        x1 = max(0, x1);
        x2 = min(maxTime, x2);

        if x2 <= x1
            continue;
        end

        patch(ax, ...
            [x1 x2 x2 x1], ...
            [0.5 0.5 nElectrodes+0.5 nElectrodes+0.5], ...
            [0.35 0.65 1.00], ...
            'FaceAlpha', 0.06, ...
            'EdgeColor', [0.10 0.30 0.80], ...
            'LineWidth', 1.3, ...
            'HandleVisibility', 'off');
    end
end
xlabel(ax,'Time (s)');
grid(ax,'on'); box(ax,'on');

% title(ax, sprintf('%s | %s | %s\nNBs(window)=%d | MeanFrag=%.2f', ...
%     methodName, cellLineName, wellName, nNB, meanFrag), ...
%     'Interpreter','none');
end

function NB = getNBstruct(NBroot, meaPlate, DIV, wellName)
NB = [];
if isfield(NBroot,meaPlate) && isfield(NBroot.(meaPlate),DIV) && ...
   isfield(NBroot.(meaPlate).(DIV),wellName)
    NB = NBroot.(meaPlate).(DIV).(wellName);
end
end

function [tFR_line, yFR_line, tFR_peaks, yFR_peaks, fragPeaks, ...
    nNB, meanFrag, startColor, endColor, lineColor, peakSpec] = ...
    extractMethodFields(NB, methodName, maxTime)

tFR_line=[]; yFR_line=[]; tFR_peaks=[]; yFR_peaks=[]; fragPeaks=[];
nNB=0; meanFrag=NaN;

if strcmpi(methodName,'Tampere')
    lineColor=[0 0 1]; startColor=[0 0.5 0]; endColor=[0.85 0 0];
    peakSpec.type='filled'; peakSpec.color='m'; peakSpec.size=22; peakSpec.marker='o';

    if ~isempty(NB) && isfield(NB,'spikeTimes')
        maskFR = NB.spikeTimes >= 0 & NB.spikeTimes <= maxTime;
        if isfield(NB,'spikeratesmooth') && ~isempty(NB.spikeratesmooth)
            tFR_line=NB.spikeTimes(maskFR); yFR_line=NB.spikeratesmooth(maskFR);
            tFR_peaks=tFR_line;             yFR_peaks=yFR_line;
        elseif isfield(NB,'spikerateplot') && ~isempty(NB.spikerateplot)
            tFR_line=NB.spikeTimes(maskFR); yFR_line=NB.spikerateplot(maskFR);
            tFR_peaks=tFR_line;             yFR_peaks=yFR_line;
        end
        if isfield(NB,'fragmentPeakTimes') && ~isempty(NB.fragmentPeakTimes)
            fragPeaks=[NB.fragmentPeakTimes{:}];
            fragPeaks=fragPeaks(fragPeaks>=0 & fragPeaks<=maxTime);
        end
        if isfield(NB,'starts'), nNB=numel(NB.starts); end
        if isfield(NB,'num_fragments') && ~isempty(NB.num_fragments)
            meanFrag=mean(NB.num_fragments);
        elseif isfield(NB,'numfragments') && ~isempty(NB.numfragments)
            meanFrag=mean(NB.numfragments);
        end
    end

else  % Doorn
    lineColor=[0 0.7 1]; startColor=[0 0.5 0]; endColor=[0.85 0 0];
    peakSpec.type='filled'; peakSpec.color='m'; peakSpec.size=22; peakSpec.marker='d';

    if ~isempty(NB)
        if isfield(NB,'spikeTimes_Doorn')
            maskFR = NB.spikeTimes_Doorn >= 0 & NB.spikeTimes_Doorn <= maxTime;
            if isfield(NB,'spikeratesmooth_Doorn') && ~isempty(NB.spikeratesmooth_Doorn)
                tFR_line=NB.spikeTimes_Doorn(maskFR); yFR_line=NB.spikeratesmooth_Doorn(maskFR);
                tFR_peaks=tFR_line; yFR_peaks=yFR_line;
            elseif isfield(NB,'spikeratesmooth_detect_Doorn') && ~isempty(NB.spikeratesmooth_detect_Doorn)
                tFR_line=NB.spikeTimes_Doorn(maskFR); yFR_line=NB.spikeratesmooth_detect_Doorn(maskFR);
                tFR_peaks=tFR_line; yFR_peaks=yFR_line;
            end
        elseif isfield(NB,'spikeTimes')
            maskFR = NB.spikeTimes >= 0 & NB.spikeTimes <= maxTime;
            if isfield(NB,'spikeratesmooth_Doorn') && ~isempty(NB.spikeratesmooth_Doorn)
                tFR_line=NB.spikeTimes(maskFR); yFR_line=NB.spikeratesmooth_Doorn(maskFR);
                tFR_peaks=tFR_line; yFR_peaks=yFR_line;
            elseif isfield(NB,'spikeratesmooth') && ~isempty(NB.spikeratesmooth)
                tFR_line=NB.spikeTimes(maskFR); yFR_line=NB.spikeratesmooth(maskFR);
                tFR_peaks=tFR_line; yFR_peaks=yFR_line;
            end
        end
        if isfield(NB,'fragmentPeakTimes') && ~isempty(NB.fragmentPeakTimes)
            fragPeaks=[NB.fragmentPeakTimes{:}];
            fragPeaks=fragPeaks(fragPeaks>=0 & fragPeaks<=maxTime);
        end
        if isfield(NB,'starts'), nNB=numel(NB.starts); end
        if isfield(NB,'num_fragments') && ~isempty(NB.num_fragments)
            meanFrag=mean(NB.num_fragments);
        elseif isfield(NB,'numfragments') && ~isempty(NB.numfragments)
            meanFrag=mean(NB.numfragments);
        end
    end
end
end

function divsSorted = sortDIVs(divNames)
divNums = zeros(numel(divNames),1);
for i = 1:numel(divNames)
    tok = regexp(divNames{i},'\d+','match','once');
    divNums(i) = str2double(tok);
end
[~,idx] = sort(divNums);
divsSorted = divNames(idx);
end

function divNum = getDIVnumber(divName)
tok = regexp(divName,'\d+','match','once');
if isempty(tok), divNum=-inf; else, divNum=str2double(tok); end
end

function electrodesSorted = sortElectrodes(electrodes)
electrodes    = cellstr(electrodes(:));
electrodeNums = zeros(numel(electrodes),1);
for i = 1:numel(electrodes)
    electrodeNums(i) = str2double(electrodes{i});
end
[~,idx] = sort(electrodeNums);
electrodesSorted = electrodes(idx);
end
