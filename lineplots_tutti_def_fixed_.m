function lineplots_tutti_def_fixed_()
% lineplots_tutti_def_fixed_  Unified NB-rate and fragments/NB line plots (Tampere + Doorn).
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Tampere_DD1C_V2,
%   NB_Tampere_DD5A_V2, NB_Tampere_TUBA_V2, NB_Tampere_C10902_V2,
%   NB_Doorn_DD1C, NB_Doorn_DD5A, NB_Doorn_TUBA, NB_Doorn_C10902
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Tampere_DD1C_V2','NB_Tampere_DD5A_V2','NB_Tampere_TUBA_V2','NB_Tampere_C10902_V2','NB_Doorn_DD1C','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% LINEPLOTS UNIFIED — NB Rate + Fragments/NB
% Consistent style across all 3 analyses:
%   Section 2 — NB Rate: Tampere + Doorn (2 separate figures)
%   Section 3 — Fragments/NB: Tampere + Doorn (2 separate figures each)
%   Section 4 — NB Rate: V1 vs V2 (1 figure with 2 subplots)
% Style rules (all plots):
%   - 1 data point = 1 well (well-level)
%   - Mean ± 95% CI (NOT SEM)
%   - IQR outlier removal (k=1.5)
%   - NB rate = nNB / recordingTime_m  [NBs/min]
%   - Line color = condition (Tampere/V1 = blue, Doorn/V2 = red-orange)
%   - Marker fill color = cell line
%   - Tampere/V1 = solid line (-o)
%   - Doorn/V2   = dashed line (--s)
%   - Shading = 95% CI in condition color
% Output:
%   - plotTvD produces 2 separate figures per call:
%       <figName>_Tampere.{png,pdf}
%       <figName>_Doorn.{png,pdf}
%   - plotV1vsV2 produces 1 figure with 2 subplots (unchanged).

%% ===== SHARED SETTINGS =====
selectedDIVs    = {'DIV53','DIV62','DIV67','DIV83','DIV88'};
DIVsnumplot     = [53 60 67 83 90];
recordingTime_s = 600;
recordingTime_m = recordingTime_s / 60;
outlierK        = 1.5;
saveFigures     = true;

saveFolder_s2 = fullfile(outRoot, 'all_outputs', 'lineplots');
saveFolder_s3 = fullfile(outRoot, 'all_outputs', 'lineplots');
saveFolder_s4 = fullfile(outRoot, 'all_outputs', 'lineplots');

for f = {saveFolder_s2, saveFolder_s3, saveFolder_s4}
    if saveFigures && ~exist(f{1},'dir'), mkdir(f{1}); end
end

%% Shared colors
colorV1 = [0.80 0.20 0.15];  % Tampere / V1 (rosso)
colorV2 = [0.15 0.40 0.75];  % Doorn / V2 (blu)

CLcolors = [0.70 0.10 0.10;   % DD1C  (rosso)
            0.95 0.55 0.10;   % DD5A  (arancio)
            0.10 0.30 0.70;   % TUBA  (blu)
            0.30 0.60 0.85];  % C10902 (azzurro)

allNames = {'DD1C','DD5A','TUBA','C10902'};

%% ===== DETECT AVAILABLE CELL LINES =====
CellLineNames = {};
for i = 1:numel(allNames)
    n  = allNames{i};
    vT = sprintf('NB_Tampere_%s', n);
    vD = sprintf('NB_Doorn_%s',   n);
    if evalin('base',sprintf('exist(''%s'',''var'')',vT)) && ...
       evalin('base',sprintf('exist(''%s'',''var'')',vD))
        CellLineNames{end+1} = n; %#ok
    else
        if ~evalin('base',sprintf('exist(''%s'',''var'')',vT))
            fprintf('[SKIP] %s not found — %s excluded from TvD.\n', vT, n);
        end
        if ~evalin('base',sprintf('exist(''%s'',''var'')',vD))
            fprintf('[SKIP] %s not found — %s excluded from TvD.\n', vD, n);
        end
    end
end

CellLineNames_V2 = {};
for i = 1:numel(allNames)
    n   = allNames{i};
    vV1 = sprintf('NB_Tampere_%s',    n);
    vV2 = sprintf('NB_Tampere_%s_V2', n);
    if evalin('base',sprintf('exist(''%s'',''var'')',vV1)) && ...
       evalin('base',sprintf('exist(''%s'',''var'')',vV2))
        CellLineNames_V2{end+1} = n; %#ok
    end
end

nCL    = numel(CellLineNames);
nCL_V2 = numel(CellLineNames_V2);
nDIV   = numel(selectedDIVs);

%% ===== COLLECT DATA =====

% --- NB Rate: Tampere vs Doorn ---
DATA_T_nb = cell(nCL, nDIV);
DATA_D_nb = cell(nCL, nDIV);

for c = 1:nCL
    n    = CellLineNames{c};
    NB_T = evalin('base', sprintf('NB_Tampere_%s', n));
    NB_D = evalin('base', sprintf('NB_Doorn_%s',   n));
    for d = 1:nDIV
        DIV = selectedDIVs{d};
        cntT = collectNBcount(NB_T, DIV, n, 'Tampere');
        cntD = collectNBcount(NB_D, DIV, n, 'Doorn');
        DATA_T_nb{c,d} = removeOutliersIQR(cntT / recordingTime_m, outlierK);
        DATA_D_nb{c,d} = removeOutliersIQR(cntD / recordingTime_m, outlierK);
    end
end

% --- NB Rate: V1 vs V2 ---
DATA_V1_nb = cell(nCL_V2, nDIV);
DATA_V2_nb = cell(nCL_V2, nDIV);

for c = 1:nCL_V2
    n    = CellLineNames_V2{c};
    NB_V1 = evalin('base', sprintf('NB_Tampere_%s',    n));
    NB_V2 = evalin('base', sprintf('NB_Tampere_%s_V2', n));
    for d = 1:nDIV
        DIV = selectedDIVs{d};
        cntV1 = collectNBcount(NB_V1, DIV, n, 'V1');
        cntV2 = collectNBcount(NB_V2, DIV, n, 'V2');
        DATA_V1_nb{c,d} = removeOutliersIQR(cntV1 / recordingTime_m, outlierK);
        DATA_V2_nb{c,d} = removeOutliersIQR(cntV2 / recordingTime_m, outlierK);
    end
end

% --- Fragments: Tampere vs Doorn ---
DATA_T_fabs = cell(nCL, nDIV);
DATA_D_fabs = cell(nCL, nDIV);
DATA_T_frat = cell(nCL, nDIV);
DATA_D_frat = cell(nCL, nDIV);

for c = 1:nCL
    n    = CellLineNames{c};
    NB_T = evalin('base', sprintf('NB_Tampere_%s', n));
    NB_D = evalin('base', sprintf('NB_Doorn_%s',   n));
    for d = 1:nDIV
        DIV = selectedDIVs{d};
        [aT,rT] = collectFragData(NB_T, DIV, n, 'Tampere');
        [aD,rD] = collectFragData(NB_D, DIV, n, 'Doorn');
        DATA_T_fabs{c,d} = removeOutliersIQR(aT, outlierK);
        DATA_D_fabs{c,d} = removeOutliersIQR(aD, outlierK);
        DATA_T_frat{c,d} = removeOutliersIQR(rT, outlierK);
        DATA_D_frat{c,d} = removeOutliersIQR(rD, outlierK);
    end
end

fprintf('\n-> All data collected.\n\n');

%% ===== SECTION 2 — NB Rate: Tampere vs Doorn =====
plotTvD(DATA_T_nb, DATA_D_nb, CellLineNames, DIVsnumplot, CLcolors, colorV1, colorV2, ...
    'NB rate (NBs/min)', ...
    sprintf('NB Rate — Tampere vs Doorn  (rec=%ds)', recordingTime_s), ...
    outlierK, saveFigures, saveFolder_s2, 'lineplot_NBrate_TampereVsDoorn');

%% ===== SECTION 3 — Fragments: Tampere vs Doorn =====
plotTvD(DATA_T_fabs, DATA_D_fabs, CellLineNames, DIVsnumplot, CLcolors, colorV1, colorV2, ...
    'Total fragments per well', ...
    'Total fragment count — Tampere vs Doorn', ...
    outlierK, saveFigures, saveFolder_s3, 'lineplot_fragments_abs_TampereVsDoorn');

plotTvD(DATA_T_frat, DATA_D_frat, CellLineNames, DIVsnumplot, CLcolors, colorV1, colorV2, ...
    'Fragments / NB', ...
    'Fragmentation ratio — Tampere vs Doorn', ...
    outlierK, saveFigures, saveFolder_s3, 'lineplot_fragments_ratio_TampereVsDoorn');

%% ===== SECTION 4 — NB Rate: V1 vs V2 =====
plotV1vsV2(DATA_V1_nb, DATA_V2_nb, CellLineNames_V2, DIVsnumplot, CLcolors, colorV1, colorV2, ...
    'NB rate (NBs/min)', ...
    sprintf('NB Rate — Tampere V1 vs V2  (rec=%ds)', recordingTime_s), ...
    outlierK, saveFigures, saveFolder_s4, 'lineplot_NBrate_V1vsV2');

fprintf('\n=== All plots saved ===\n');

%% PLOT FUNCTIONS

end

function plotTvD(DATA_T, DATA_D, CellLineNames, DIVs, CLcolors, colorT, colorD, ...
    yLabelStr, titleStr, outlierK, saveFig, saveFolder, figName)
% Produces TWO separate figures:
%   1) Tampere only — saved as [figName '_Tampere']
%   2) Doorn   only — saved as [figName '_Doorn']
% Each figure shows all cell lines on the same axes:
%   - Tampere figure: solid line (-o), line color = colorT
%   - Doorn   figure: dashed line (--s), line color = colorD
%   - Marker fill color = cell line (CLcolors)
%   - Shaded band = 95% CI in condition color

nCL = numel(CellLineNames); nDIV = numel(DIVs);
[muT,ciT] = computeMeanCI(DATA_T, nCL, nDIV);
[muD,ciD] = computeMeanCI(DATA_D, nCL, nDIV);

% --- Figure 1: Tampere only ---
plotSingleDetector(muT, ciT, CellLineNames, DIVs, CLcolors, colorT, ...
    '-o', yLabelStr, [titleStr ' — Tampere'], outlierK, ...
    saveFig, saveFolder, [figName '_Tampere']);

% --- Figure 2: Doorn only ---
plotSingleDetector(muD, ciD, CellLineNames, DIVs, CLcolors, colorD, ...
    '--s', yLabelStr, [titleStr ' — Doorn'], outlierK, ...
    saveFig, saveFolder, [figName '_Doorn']);
end

function plotSingleDetector(mu, ci, CellLineNames, DIVs, CLcolors, lineColor, ...
    lineStyle, yLabelStr, titleStr, outlierK, saveFig, saveFolder, figName)
% Helper: draw one figure for one detector across all cell lines.
% Style (matches script2_lineplots_new):
%   - Line color    = cell line (CLcolors)
%   - Marker        = circle (-o) for all cell lines, uniform
%   - Marker fill   = same as line color (filled circle)
%   - CI            = SHADED band (no error bars)
%   - Grid          = major only
%   - Title         = "<title> | IQR k=<k>"
%   - Ylabel includes "± 95% CI"
% Note: lineColor and lineStyle are kept in the signature for
% backward compatibility with plotTvD but are NOT used here -- each
% cell line uses its own color from CLcolors, and the marker is
% always '-o'. The Tampere / Doorn distinction is conveyed by the
% figure title and file name, not by linestyle.

nCL = numel(CellLineNames);
fig = figure('Color','w','Position',[100 100 800 520]);
ax  = gca;
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
hLines = gobjects(1,nCL);

for c = 1:nCL
    col = CLcolors(c,:);
    [dv,mv,ev] = getValidPoints(DIVs, mu(c,:), ci(c,:));
    if ~isempty(dv)
        % Shaded CI band
        xfill = [dv, fliplr(dv)];
        yfill = [mv+ev, fliplr(mv-ev)];
        fill(ax, xfill, yfill, col, 'FaceAlpha', 0.15, ...
             'EdgeColor', 'none', 'HandleVisibility', 'off');
        % Line + filled circle markers
        hLines(c) = plot(ax, dv, mv, '-o', ...
            'Color',           col, ...
            'LineWidth',       1.8, ...
            'MarkerSize',      6, ...
            'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col);
    else
        hLines(c) = plot(ax, nan, nan, '-o', ...
            'Color', col, 'MarkerFaceColor', col);
    end
end

% --- Axis style (image-matching version) ---
set(ax, 'XTick', DIVs, ...
    'XTickLabel', arrayfun(@(x) sprintf('DIV%d',x), DIVs, 'UniformOutput', false), ...
    'FontSize', 10, 'LineWidth', 0.8);
xlabel(ax, 'DIV', 'FontSize', 15);
ylabel(ax, [yLabelStr ' \pm 95% CI'], 'FontSize', 15);
% title(ax, sprintf('%s | IQR outlier removal (k=%.1f)', titleStr, outlierK), ...
%     'FontSize', 12, 'FontWeight', 'bold');
xlim(ax, [DIVs(1)-4, DIVs(end)+4]);

% --- Legend with box (matches image) ---
legend(hLines, CellLineNames, 'Location','northwest', ...
    'FontSize',10, 'Box','on');

saveAndClose(fig, saveFig, saveFolder, figName);
end

function plotV1vsV2(DATA_V1, DATA_V2, CellLineNames, DIVs, CLcolors, colorV1, colorV2, ...
    yLabelStr, titleStr, outlierK, saveFig, saveFolder, figName)
% 2 subplots side by side: left = V1, right = V2
% Each subplot: 1 line per cell line, cell line color

nCL = numel(CellLineNames); nDIV = numel(DIVs);
[muV1,ciV1] = computeMeanCI(DATA_V1, nCL, nDIV);
[muV2,ciV2] = computeMeanCI(DATA_V2, nCL, nDIV);

fig = figure('Color','w','Position',[60 60 1200 520]);
tl  = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% title(tl, sprintf('%s  (IQR k=%.1f | mean ± 95%% CI)', titleStr, outlierK), ...
%     'FontSize',12,'FontWeight','bold');

datasets  = {muV1, ciV1, muV2, ciV2};
colors_   = {colorV1, colorV2};
subTitles = {'Tampere V1', 'Tampere V2'};

for v = 1:2
    ax  = nexttile;
    mu_ = datasets{(v-1)*2+1};
    ci_ = datasets{(v-1)*2+2};
    col_ = colors_{v};

    hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    hLines = gobjects(1,nCL);

    for c = 1:nCL
        col = CLcolors(c,:);
        [dv,mv,ev] = getValidPoints(DIVs, mu_(c,:), ci_(c,:));
        if ~isempty(dv)
            fill(ax,[dv fliplr(dv)],[mv+ev fliplr(max(mv-ev,0))], col, ...
                'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
            hTemp = plot(ax,dv,mv,'-o','Color',col,'LineWidth',1.8, ...
                'MarkerSize',6,'MarkerFaceColor',col,'MarkerEdgeColor',col);
            hLines(c) = hTemp(1);
        else
            hLines(c) = plot(ax,nan,nan,'-o','Color',col,'MarkerFaceColor',col);
        end
    end

    applyAxisStyle(ax, DIVs, yLabelStr, subTitles{v}, outlierK);
    legend(ax, hLines, CellLineNames, 'Location','northwest','FontSize',10,'Box','off');
end

saveAndClose(fig, saveFig, saveFolder, figName);
end

%% DATA COLLECTION FUNCTIONS

function vals = collectNBcount(NB_Cell, DIV, clName, detName)
vals = [];
if isempty(NB_Cell), return; end
plateNames = fieldnames(NB_Cell);
for p = 1:numel(plateNames)
    meaPlate = plateNames{p};
    if ~isfield(NB_Cell.(meaPlate),DIV), continue; end
    wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));
    for w = 1:numel(wellNames)
        NB = NB_Cell.(meaPlate).(DIV).(wellNames{w});
        nNB = 0;
        if isstruct(NB) && isfield(NB,'starts') && ~isempty(NB.starts) && ~all(isnan(NB.starts))
            nNB = numel(NB.starts);
        elseif isstruct(NB) && isfield(NB,'num_bursts') && ~isempty(NB.num_bursts)
            nNB = NB.num_bursts;
        elseif isnumeric(NB) && size(NB,2) >= 2
            nNB = size(NB,1);
        end
        if nNB == 0
            fprintf('[INACTIVE] %s | %s | %s | %s — excluded.\n', ...
                clName, detName, DIV, wellNames{w});
            continue
        end
        vals(end+1) = nNB; %#ok
    end
end
end

function [absVals, ratVals] = collectFragData(NB_Cell, DIV, clName, detName)
absVals = []; ratVals = [];
if isempty(NB_Cell) || ~isstruct(NB_Cell), return; end
plateNames = fieldnames(NB_Cell);
for p = 1:numel(plateNames)
    meaPlate = plateNames{p};
    if ~isfield(NB_Cell.(meaPlate),DIV), continue; end
    wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));
    for w = 1:numel(wellNames)
        NB = NB_Cell.(meaPlate).(DIV).(wellNames{w});
        if ~isstruct(NB) || ~isfield(NB,'starts') || ...
           isempty(NB.starts) || all(isnan(NB.starts)), continue; end
        if ~isfield(NB,'num_fragments') || isempty(NB.num_fragments), continue; end
        frags = NB.num_fragments(:); frags = frags(isfinite(frags));
        frags_valid = frags(frags >= 1);
        if isempty(frags_valid), continue; end
        absVals(end+1) = sum(frags_valid);                      %#ok
        ratVals(end+1) = sum(frags_valid) / numel(frags_valid); %#ok
    end
end
end

%% UTILITY FUNCTIONS

function [mu, ci] = computeMeanCI(DATA, nCL, nDIV)
mu  = nan(nCL,nDIV);
ci  = nan(nCL,nDIV);
for c = 1:nCL
    for d = 1:nDIV
        v = DATA{c,d};
        if numel(v) >= 2
            mu(c,d) = mean(v);
            ci(c,d) = tinv(0.975, numel(v)-1) * std(v) / sqrt(numel(v));
        end
    end
end
end

function [dv, mv, ev] = getValidPoints(DIVs, mu_row, ci_row)
vi = ~isnan(mu_row);
dv = DIVs(vi); mv = mu_row(vi); ev = ci_row(vi);
end

function applyAxisStyle(ax, DIVs, yLabelStr, titleStr, outlierK)
set(ax, 'XTick', DIVs, ...
    'XTickLabel', arrayfun(@(x) sprintf('DIV%d',x), DIVs, 'UniformOutput', false), ...
    'FontSize', 10, 'LineWidth', 0.8);
xlabel(ax, 'DIV', 'FontSize', 15, 'FontWeight', 'bold');
ylabel(ax, yLabelStr, 'FontSize', 15);
% title(ax, titleStr, 'FontSize', 12, 'FontWeight', 'bold');
subtitle(ax, sprintf('IQR k=%.1f | mean ± 95%% CI | 1 point = 1 well', outlierK), ...
    'FontSize', 9, 'Color', [0.4 0.4 0.4]);
xlim(ax, [DIVs(1)-4, DIVs(end)+4]);
end

function saveAndClose(fig, saveFig, saveFolder, figName)
if saveFig
    exportgraphics(fig, fullfile(saveFolder,[figName '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(saveFolder,[figName '.pdf']), 'ContentType', 'vector');
    fprintf('Saved: %s\n', figName);
end
close(fig);
end

function clean = removeOutliersIQR(data, k)
data = data(isfinite(data));
if numel(data) < 4, clean = data; return; end
q1   = prctile(data, 25);
q3   = prctile(data, 75);
iqr_ = q3 - q1;
clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end
