function RASTER_NBandSEB()
% RASTER_NBandSEB  Batch raster figures for single-electrode bursts and network bursts (Doorn/Tampere).
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Doorn_DD1C, SampFreq,
%   ParamsBurst, maxTime, timeBin, movingWindow_sec, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Doorn_DD5A, NB_Doorn_TUBA,
%   NB_Doorn_C10902
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Doorn_DD1C','SampFreq','ParamsBurst','maxTime','timeBin','movingWindow_sec','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

%% NO RICALCOLO NB
% BATCH RASTER SCRIPT
% Creates one figure for each cell line × MEA plate
% NB Tampere: read from NB_Tampere_<CellLine> (pre-computed)
% NB Doorn:   read from NB_Doorn_<CellLine>   (pre-computed)
% SEB:        ricalcola al volo da spikeTable + DetectBurstsNeuro
% REQUIRED IN WORKSPACE:
%   DD1C, DD5A, TUBA, C10902      (raw spike data)
%   NB_Tampere_DD1C, ...          (pre-computed Tampere NB structures)
%   NB_Doorn_DD1C, ...            (pre-computed Doorn NB structures)
%   ParamsBurst

%% ===== COMMON CELL LINE LIST =====
CellLines = buildCellLines();
if isempty(CellLines)
    error('No valid cell line variables found in workspace.');
end

% 1) SINGLE-ELECTRODE BURST RASTER - 4x2 grid (rows = cell lines, cols = DIVs)
cfgSE             = struct();
cfgSE.cellOrder   = {'DD1C','DD5A','TUBA','C10902'};   % rows
cfgSE.DIVs        = {'DIV53','DIV90'};                 % columns
cfgSE.maxTime     = 100;                               % s shown per panel
cfgSE.mergeGapSec = 0.13;
cfgSE.SampFreq    = 12500;
cfgSE.recDuration = 600;
cfgSE.saveFigure  = true;
cfgSE.saveFolder  = fullfile(outRoot, 'all_outputs', 'raster_SEB');
% Optional explicit well per cell line, e.g.
%   cfgSE.wellPick.DD1C = struct('plate','N2894','well','B2');
% Otherwise the first well valid at BOTH DIVs is auto-selected.
plotSingleElectrodeGrid(CellLines, cfgSE);

% 2) NETWORK BURST RASTER GRID — rows = Doorn/Tampere, cols = cell line x DIV
%    Two figures: one for the DS lines, one for the control lines.
cfgNB             = struct();
cfgNB.DIVs        = {'DIV53','DIV90'};                  % columns (per cell line)
cfgNB.maxTime     = 120;                                % s shown per panel
cfgNB.recDuration = 600;
cfgNB.SampFreq    = 12500;
cfgNB.saveFigure  = true;
cfgNB.saveFolder  = fullfile(outRoot, 'all_outputs', 'raster_NB');
cfgNB.ParamsDoorn = struct('timeBin',25e-3);            % only for numeric-format Doorn NB
cfgNB.groups      = struct('DS',{{'DD1C','DD5A'}}, 'Control',{{'TUBA','C10902'}});
% Optional explicit well per cell line:
%   cfgNB.wellPick.DD1C = struct('plate','N2894','well','B2');

CellLines_Doorn   = buildCellLinesNB('Doorn');
CellLines_Tampere = buildCellLinesNB('Tampere');
if isempty(CellLines_Doorn) && isempty(CellLines_Tampere)
    warning('No NB_Doorn_* or NB_Tampere_* structures found — NB raster grid skipped.');
else
    plotNBrasterGrid(CellLines, CellLines_Tampere, CellLines_Doorn, cfgNB);
end

% LOCAL FUNCTIONS

% Build cell-line list with raw spike data
end

function CellLines = buildCellLines()
CellLines = {};
names = {'DD1C','DD5A','TUBA','C10902'};
for i = 1:numel(names)
    n = names{i};
    if evalin('base', sprintf('exist(''%s'',''var'')', n))
        CellLines{end+1} = struct('name', n, 'Data', evalin('base', n)); %#ok<AGROW>
    end
end
end

% Build cell-line list with the pre-computed NB structures
function CellLines_NB = buildCellLinesNB(detector)
% detector: 'Tampere' oppure 'Doorn'
CellLines_NB = {};
names = {'DD1C','DD5A','TUBA','C10902'};
for i = 1:numel(names)
    n   = names{i};
    var = sprintf('NB_%s_%s', detector, n);
    if evalin('base', sprintf('exist(''%s'',''var'')', var))
        CellLines_NB{end+1} = struct( ...
            'name', n, ...
            'NB',   evalin('base', var)); %#ok<AGROW>
    end
end
end

% SINGLE-ELECTRODE BURSTS

function plotSingleElectrodeGrid(CellLines, cfg)
% plotSingleElectrodeGrid  One figure, 4 x 2 single-electrode burst rasters:
%   rows    = cfg.cellOrder (cell lines), labelled on the left
%   columns = cfg.DIVs,      labelled on top
% One representative well per cell line is shown at both DIVs (same well).
if cfg.saveFigure && ~exist(cfg.saveFolder,'dir'), mkdir(cfg.saveFolder); end

% Map cell line name -> raw data struct
nameMap = containers.Map('KeyType','char','ValueType','any');
for c = 1:numel(CellLines)
    nameMap(CellLines{c}.name) = CellLines{c}.Data;
end

cellOrder = cfg.cellOrder;
DIVs      = cfg.DIVs;
nR = numel(cellOrder);
nC = numel(DIVs);

fig = figure('Color','w','Name','Single-electrode burst rasters (grid)', ...
    'NumberTitle','off','Position',[40 40 560*nC 330*nR]);
tl = tiledlayout(nR, nC, 'TileSpacing','compact','Padding','compact');

for r = 1:nR
    clName = cellOrder{r};
    if ~isKey(nameMap, clName)
        for k = 1:nC, ax = nexttile; axis(ax,'off'); end
        continue
    end
    CellData = nameMap(clName);
    [plateSel, wellSel] = pickWellAcrossDIVs(CellData, DIVs, cfg, clName);

    for k = 1:nC
        ax = nexttile; hold(ax,'on');
        DIV = DIVs{k};
        if ~isempty(plateSel) && isfield(CellData,plateSel) && ...
           isfield(CellData.(plateSel),DIV) && isfield(CellData.(plateSel).(DIV),wellSel)
            drawSingleElectrodeRaster(ax, CellData.(plateSel).(DIV).(wellSel), cfg);
        else
            axis(ax,'off');
        end
        if r == 1   % column header = DIV
%             title(ax, strrep(DIV,'DIV','DIV '), 'FontSize', 16, 'FontWeight','bold');
        end
        if k == 1   % row label = cell line
            ylabel(ax, clName, 'FontSize', 15, 'FontWeight','bold', 'Interpreter','none');
        else
            ylabel(ax, '');
        end
        if r == nR
            xlabel(ax, 'Time (s)', 'FontSize', 15);
        else
            xlabel(ax, '');
        end
    end
end

% sgtitle(tl, 'Single-electrode burst rasters', 'FontSize', 16, 'FontWeight','bold');

if cfg.saveFigure
    exportgraphics(fig, fullfile(cfg.saveFolder, 'singleElectrodeRaster_grid_4x2.png'), ...
        'Resolution', 300);
end
end

function [plateSel, wellSel] = pickWellAcrossDIVs(CellData, DIVs, cfg, clName)
% Choose one well (and its plate) that has valid spikes at ALL requested DIVs.
% An explicit choice in cfg.wellPick.(clName) = struct('plate',..,'well',..) wins.
plateSel = ''; wellSel = '';
if isfield(cfg,'wellPick') && isfield(cfg.wellPick, clName)
    plateSel = cfg.wellPick.(clName).plate;
    wellSel  = cfg.wellPick.(clName).well;
    return
end
plates = fieldnames(CellData);
for p = 1:numel(plates)
    pl = plates{p};
    if ~all(cellfun(@(d) isfield(CellData.(pl),d), DIVs)), continue; end
    wsets  = cellfun(@(d) fieldnames(CellData.(pl).(d)), DIVs, 'UniformOutput', false);
    common = wsets{1};
    for q = 2:numel(wsets), common = intersect(common, wsets{q}); end
    for wi = 1:numel(common)
        wn = common{wi};
        ok = true;
        for q = 1:numel(DIVs)
            st = CellData.(pl).(DIVs{q}).(wn);
            if ~istable(st) || height(st)==0 || ~ismember('Time',st.Properties.VariableNames)
                ok = false; break
            end
        end
        if ok, plateSel = pl; wellSel = wn; return; end
    end
end
end

function drawSingleElectrodeRaster(ax, spikeTable, cfg)
% Draw one single-electrode raster (spikes as ticks, detected bursts as red
% patches) for a single well into axes ax.
if ~istable(spikeTable) || height(spikeTable)==0 || ...
   ~ismember('Channel',spikeTable.Properties.VariableNames) || ...
   ~ismember('Time',   spikeTable.Properties.VariableNames)
    axis(ax,'off'); return
end
Time    = spikeTable.Time;
Channel = string(spikeTable.Channel);
valid   = isfinite(Time) & Time >= 0;
Time    = Time(valid); Channel = Channel(valid);
if isempty(Time), axis(ax,'off'); return; end

t0       = min(Time);
Time0    = Time - t0;
showEnd  = min(cfg.maxTime, cfg.recDuration);
maskShow = Time0 >= 0 & Time0 <= showEnd;
TimeShow    = Time0(maskShow);
ChannelShow = Channel(maskShow);
if isempty(TimeShow), axis(ax,'off'); return; end

electrodes  = sortElectrodesSafe(unique(ChannelShow,'stable'));
nElectrodes = numel(electrodes);

for e = 1:nElectrodes
    elec        = string(electrodes{e});
    elecSpkShow = TimeShow(ChannelShow == elec);
    plot(ax, elecSpkShow, e*ones(numel(elecSpkShow),1), '|', ...
        'Color',[0 0 0],'MarkerSize',6);

    spT = sort(Time0(Channel == elec));
    sps = unique(round(spT * cfg.SampFreq) + 1);
    sps = sps(sps >= 1);
    if numel(sps) < 3, continue; end

    bt = DetectBurstsNeuro(sps, cfg.SampFreq, cfg.recDuration, evalin('base','ParamsBurst'));
    if isempty(bt), continue; end
    T_s = (bt(:,1)-1)/cfg.SampFreq;
    T_e = (bt(:,2)-1)/cfg.SampFreq;
    if size(bt,1) > 1
        [T_s, T_e] = MergeBurstsUniversal([], T_s, T_e, cfg.mergeGapSec, 'Static');
    end
    inW  = T_s <= showEnd & T_e >= 0;
    bSec = T_s(inW); eSec = T_e(inW);
    for b = 1:numel(bSec)
        x1 = max(0, bSec(b));
        x2 = min(showEnd, eSec(b));
        if x2 <= x1, continue; end
        patch(ax, [x1 x2 x2 x1], [e-0.38 e-0.38 e+0.38 e+0.38], [1.00 0.35 0.35], ...
            'FaceAlpha',0.10,'EdgeColor',[0.80 0.10 0.10],'LineWidth',1.0, ...
            'HandleVisibility','off');
    end
end

ylim(ax,[0 nElectrodes+1]); xlim(ax,[0 showEnd]);
yticks(ax,1:nElectrodes); yticklabels(ax,electrodes);
set(ax,'FontSize',14); grid(ax,'on'); box(ax,'on');
end

% NETWORK BURST RASTERS — grid (Doorn + Tampere)

function plotNBrasterGrid(CellLines, CellLines_Tampere, CellLines_Doorn, cfg)
% plotNBrasterGrid  Two NB-raster figures (one per cell-line group). Each is a
%   2 x (2*nLines) grid: rows = Doorn (top) / Tampere (bottom), columns =
%   cell line x DIV. One representative well per cell line is shown.
if cfg.saveFigure && ~exist(cfg.saveFolder,'dir'), mkdir(cfg.saveFolder); end

dataMap = mapByName(CellLines);
tMap    = mapByName(CellLines_Tampere);
dMap    = mapByName(CellLines_Doorn);

DIVs       = cfg.DIVs;
detectors  = {'Doorn','Tampere'};
groupNames = fieldnames(cfg.groups);

for g = 1:numel(groupNames)
    gname = groupNames{g};
    lines = cfg.groups.(gname);
    nCL   = numel(lines);
    nC    = nCL * numel(DIVs);

    % One representative well per cell line (valid spikes + NB present at both DIVs)
    pickPlate = cell(1,nCL); pickWell = cell(1,nCL);
    for ci = 1:nCL
        [pickPlate{ci}, pickWell{ci}] = pickWellForNB( ...
            getMapVal(dataMap,lines{ci}), getMapVal(tMap,lines{ci}), ...
            getMapVal(dMap,lines{ci}), DIVs, cfg, lines{ci});
    end

    fig = figure('Color','w','Name',sprintf('NB rasters - %s lines',gname), ...
        'NumberTitle','off','Position',[30 30 max(900,460*nC) 720]);
    tl = tiledlayout(2, nC, 'TileSpacing','compact','Padding','compact');

    for dr = 1:numel(detectors)
        det = detectors{dr};
        for ci = 1:nCL
            clName   = lines{ci};
            CellData = getMapVal(dataMap, clName);
            if strcmp(det,'Doorn'), NBroot = getMapVal(dMap, clName);
            else,                   NBroot = getMapVal(tMap, clName); end
            pl = pickPlate{ci}; wl = pickWell{ci};

            for k = 1:numel(DIVs)
                DIV = DIVs{k};
                col = (ci-1)*numel(DIVs) + k;
                ax  = nexttile((dr-1)*nC + col); hold(ax,'on');

                ok = ~isempty(pl) && ~isempty(CellData) && isfield(CellData,pl) && ...
                     isfield(CellData.(pl),DIV) && isfield(CellData.(pl).(DIV),wl);
                if ok
                    spikeTable = CellData.(pl).(DIV).(wl);
                    NBwell = [];
                    if ~isempty(NBroot) && isfield(NBroot,pl) && isfield(NBroot.(pl),DIV) ...
                            && isfield(NBroot.(pl).(DIV),wl)
                        NBwell = NBroot.(pl).(DIV).(wl);
                    end
                    drawNBraster(ax, spikeTable, NBwell, det, cfg);
                else
                    axis(ax,'off');
                end

                if dr == 1
%                     title(ax, sprintf('%s - DIV %s', clName, erase(DIV,'DIV')), ...
%                         'FontSize',14,'FontWeight','bold','Interpreter','none');
                end
                if col == 1
                    ylabel(ax, det, 'FontSize',15,'FontWeight','bold');
                else
                    ylabel(ax, '');
                end
                if dr == numel(detectors)
                    xlabel(ax, 'Time (s)', 'FontSize',15);
                else
                    xlabel(ax, '');
                end
            end
        end
    end

%     sgtitle(tl, sprintf('Network-burst rasters - %s lines (first %d s)', gname, cfg.maxTime), ...
%         'FontSize',15,'FontWeight','bold');
    if cfg.saveFigure
        exportgraphics(fig, fullfile(cfg.saveFolder, sprintf('NBraster_grid_%s.png',gname)), ...
            'Resolution',300);
    end
end
end

function drawNBraster(ax, spikeTable, NBwell, detector, cfg)
% Draw one NB raster (spike ticks + NB patches) for a single well/detector.
if ~istable(spikeTable) || height(spikeTable)==0 || ...
   ~ismember('Channel',spikeTable.Properties.VariableNames) || ...
   ~ismember('Time',   spikeTable.Properties.VariableNames)
    axis(ax,'off'); return
end
Time    = spikeTable.Time;
Channel = string(spikeTable.Channel);
valid   = isfinite(Time) & Time >= 0;
Time    = Time(valid); Channel = Channel(valid);
if isempty(Time), axis(ax,'off'); return; end

t0       = min(Time);
Time0    = Time - t0;
showEnd  = min(cfg.maxTime, cfg.recDuration);
maskShow = Time0 >= 0 & Time0 <= showEnd;
TimeShow    = Time0(maskShow);
ChannelShow = Channel(maskShow);

electrodes  = sortElectrodesSafe(unique(ChannelShow,'stable'));
nElectrodes = numel(electrodes);
for e = 1:nElectrodes
    elec        = string(electrodes{e});
    elecSpkShow = TimeShow(ChannelShow == elec);
    plot(ax, elecSpkShow, e*ones(numel(elecSpkShow),1), '|', 'Color',[0 0 0],'MarkerSize',6);
end

if strcmp(detector,'Doorn')
    faceC = [0.35 0.65 1.00]; edgeC = [0.10 0.30 0.80];
else
    faceC = [1.00 0.45 0.45]; edgeC = [0.80 0.10 0.10];
end

[nbStarts, nbEnds] = extractNBwindows(NBwell, cfg);
nbShown = 0;
for b = 1:numel(nbStarts)
    x1 = max(0, nbStarts(b));
    x2 = min(showEnd, nbEnds(b));
    if x2 <= x1, continue; end
    patch(ax, [x1 x2 x2 x1], [0.5 0.5 nElectrodes+0.5 nElectrodes+0.5], faceC, ...
        'FaceAlpha',0.07,'EdgeColor',edgeC,'LineWidth',1.2,'HandleVisibility','off');
    nbShown = nbShown + 1;
end

ylim(ax,[0 nElectrodes+1]); xlim(ax,[0 showEnd]);
yticks(ax,1:nElectrodes); yticklabels(ax,electrodes);
set(ax,'FontSize',14); grid(ax,'on'); box(ax,'on');
text(ax, showEnd, nElectrodes+0.7, sprintf('NBs=%d',nbShown), ...
    'HorizontalAlignment','right','VerticalAlignment','top','FontSize',9,'Color',edgeC);
end

function [nbStarts, nbEnds] = extractNBwindows(NB, cfg)
% NB start/end times [s] from a Tampere/Doorn NB struct or numeric matrix.
nbStarts = []; nbEnds = [];
if isempty(NB), return; end
if isstruct(NB) && isfield(NB,'starts') && isfield(NB,'ends') && ~isempty(NB.starts)
    nbStarts = NB.starts(:);
    nbEnds   = NB.ends(:);
elseif isnumeric(NB) && size(NB,2) >= 2
    tb = 25e-3;
    if isfield(cfg,'ParamsDoorn') && isfield(cfg.ParamsDoorn,'timeBin'), tb = cfg.ParamsDoorn.timeBin; end
    nbStarts = (NB(:,1)-1) * tb;
    nbEnds   = (NB(:,2)-1) * tb;
end
end

function [plateSel, wellSel] = pickWellForNB(CellData, NB_T, NB_D, DIVs, cfg, clName)
% Pick a well with valid spikes at all DIVs and an NB entry in at least one
% detector. Explicit cfg.wellPick.(clName) wins; falls back to any valid well.
plateSel = ''; wellSel = '';
if isfield(cfg,'wellPick') && isfield(cfg.wellPick, clName)
    plateSel = cfg.wellPick.(clName).plate;
    wellSel  = cfg.wellPick.(clName).well;
    return
end
if isempty(CellData), return; end
plates = fieldnames(CellData);
for p = 1:numel(plates)
    pl = plates{p};
    if ~all(cellfun(@(d) isfield(CellData.(pl),d), DIVs)), continue; end
    wsets  = cellfun(@(d) fieldnames(CellData.(pl).(d)), DIVs, 'UniformOutput', false);
    common = wsets{1};
    for q = 2:numel(wsets), common = intersect(common, wsets{q}); end
    for wi = 1:numel(common)
        wn = common{wi};
        okData = true;
        for q = 1:numel(DIVs)
            st = CellData.(pl).(DIVs{q}).(wn);
            if ~istable(st) || height(st)==0 || ~ismember('Time',st.Properties.VariableNames)
                okData = false; break
            end
        end
        if ~okData, continue; end
        if nbHasWell(NB_T, pl, DIVs, wn) || nbHasWell(NB_D, pl, DIVs, wn)
            plateSel = pl; wellSel = wn; return
        end
    end
end
% fallback: first data-valid well even without an NB entry
for p = 1:numel(plates)
    pl = plates{p};
    if ~all(cellfun(@(d) isfield(CellData.(pl),d), DIVs)), continue; end
    wsets  = cellfun(@(d) fieldnames(CellData.(pl).(d)), DIVs, 'UniformOutput', false);
    common = wsets{1};
    for q = 2:numel(wsets), common = intersect(common, wsets{q}); end
    if ~isempty(common), plateSel = pl; wellSel = common{1}; return; end
end
end

function tf = nbHasWell(NBroot, plate, DIVs, well)
tf = false;
if isempty(NBroot) || ~isfield(NBroot,plate), return; end
for q = 1:numel(DIVs)
    if ~isfield(NBroot.(plate),DIVs{q}) || ~isfield(NBroot.(plate).(DIVs{q}),well)
        return
    end
end
tf = true;
end

function m = mapByName(CellLinesCell)
% name -> payload map from a cell array of structs with .name and (.Data | .NB)
m = struct('keys',{{}},'vals',{{}});
for i = 1:numel(CellLinesCell)
    sline = CellLinesCell{i};
    m.keys{end+1} = sline.name;
    if isfield(sline,'Data'),   m.vals{end+1} = sline.Data;
    elseif isfield(sline,'NB'), m.vals{end+1} = sline.NB;
    else,                       m.vals{end+1} = []; end
end
end

function v = getMapVal(m, name)
v = [];
idx = find(strcmp(m.keys, name), 1);
if ~isempty(idx), v = m.vals{idx}; end
end

function electrodesSorted = sortElectrodesSafe(electrodes)
electrodes = cellstr(string(electrodes(:)));
numVal     = nan(numel(electrodes),1);
for i = 1:numel(electrodes)
    tok = regexp(electrodes{i},'\d+','match','once');
    if ~isempty(tok), numVal(i) = str2double(tok); end
end
[~, idx]        = sortrows([isnan(numVal) numVal],[1 2]);
electrodesSorted = electrodes(idx);
end
