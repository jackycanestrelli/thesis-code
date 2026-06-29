function boxplot_TvD_fragcomaprison()
% boxplot_TvD_fragcomaprison  Fragments-per-NB boxplots (absolute and ratio) comparing Tampere vs Doorn.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Doorn_DD1C, NB_Doorn_DD5A,
%   NB_Doorn_TUBA, NB_Doorn_C10902
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Doorn_DD1C','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% BOXPLOT: TAMPERE vs DOORN — Fragments comparison per cell line
% Data point:
%   Fig 1 (abs)   = total num of fragments per well (sum over all NBs)
%                   → 1 value per well
%   Fig 2 (ratio) = sum(numfragments_valid) / num_NB_with_frags per well
%                   → 1 value per well
%                   NB with numfragments == 0 are EXCLUDED from both
%                   numerator and denominator (debug printed).
% REQUIRED IN WORKSPACE:
%   NB_Tampere_DD1C, NB_Tampere_DD5A, NB_Tampere_TUBA, NB_Tampere_C10902
%   NB_Doorn_DD1C,   NB_Doorn_DD5A,   NB_Doorn_TUBA,   NB_Doorn_C10902

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
outlierK     = 1.5;
saveFigures  = true;
saveFolder   = fullfile(outRoot, 'all_outputs', 'frag_TvsD');

if saveFigures && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

colorT = [0.80 0.20 0.15];
colorD = [0.15 0.40 0.75];

%% ===== CELL LINE DEFINITIONS =====
allNames      = {'DD1C','DD5A','TUBA','C10902'};
CellLineNames = {};
for i = 1:numel(allNames)
    n  = allNames{i};
    vT = sprintf('NB_Tampere_%s', n);
    vD = sprintf('NB_Doorn_%s',   n);
    if exist(vT,'var') && exist(vD,'var')
        CellLineNames{end+1} = n; %#ok
    else
        if ~exist(vT,'var'), fprintf('[SKIP] %s not found — %s excluded.\n', vT, n); end
        if ~exist(vD,'var'), fprintf('[SKIP] %s not found — %s excluded.\n', vD, n); end
    end
end
if isempty(CellLineNames)
    error('No valid cell line with both NB_Tampere_* and NB_Doorn_* found.');
end

nCL  = numel(CellLineNames);
nDIV = numel(selectedDIVs);

%% ===== DATA COLLECTION =====
DATA_T_abs = cell(nCL, nDIV);
DATA_D_abs = cell(nCL, nDIV);
DATA_T_rat = cell(nCL, nDIV);
DATA_D_rat = cell(nCL, nDIV);

for c = 1:nCL
    n    = CellLineNames{c};
    NB_T = eval(sprintf('NB_Tampere_%s', n));
    NB_D = eval(sprintf('NB_Doorn_%s',   n));

    for d = 1:nDIV
        DIV = selectedDIVs{d};

        [absT, ratT] = collectFragData(NB_T, DIV, n, 'Tampere', DIV);
        [absD, ratD] = collectFragData(NB_D, DIV, n, 'Doorn',   DIV);

        if outlierK > 0
            absT = removeOutliersIQR(absT, outlierK);
            absD = removeOutliersIQR(absD, outlierK);
            ratT = removeOutliersIQR(ratT, outlierK);
            ratD = removeOutliersIQR(ratD, outlierK);
        end

        DATA_T_abs{c,d} = absT;
        DATA_D_abs{c,d} = absD;
        DATA_T_rat{c,d} = ratT;
        DATA_D_rat{c,d} = ratD;

        fprintf('[DATA] %s | %s | Tampere: %d wells (abs), %d wells (ratio) | Doorn: %d wells (abs), %d wells (ratio)\n', ...
            n, DIV, numel(absT), numel(ratT), numel(absD), numel(ratD));
    end
end

%% ===== FIGURE 1 — Total fragments per well =====
plotDivGrouped(DATA_T_abs, DATA_D_abs, CellLineNames, selectedDIVs, ...
    colorT, colorD, outlierK, ...
    'Total number of fragments (per well)', ...
    'Total fragment count per well — Tampere vs Doorn', ...
    'boxplot_fragments_abs_TampereVsDoorn', ...
    '1 data point = 1 well  |  sum of numfragments (NB with frags >= 1 only)', ...
    saveFigures, saveFolder);

%% ===== FIGURE 2 — Fragments/NB ratio per well =====
plotDivGrouped(DATA_T_rat, DATA_D_rat, CellLineNames, selectedDIVs, ...
    colorT, colorD, outlierK, ...
    'Fragments / NB  (NBs/min)', ...
    'Fragmentation ratio — Tampere vs Doorn | NBs with at least 1 frag', ...
    'boxplot_fragments_ratio_TampereVsDoorn', ...
    '1 data point = 1 well  ', ...
    saveFigures, saveFolder);

% LOCAL FUNCTIONS

end

function [absVals, ratVals] = collectFragData(NB_Cell, DIV, clName, detName, divName)
% Per ogni well attivo:
%   absVals: sum(numfragments) per i soli NB con numfragments >= 1
%   ratVals: sum(numfragments_valid) / numel(NB_valid) — ratio per well
% Well inattivi (starts vuoto/NaN) → [INACTIVE] debug + esclusi.
% NB con numfragments == 0          → [FRAG_SKIP] debug + esclusi dal calcolo.

    absVals = [];
    ratVals = [];
    if isempty(NB_Cell) || ~isstruct(NB_Cell), return; end

    plateNames = fieldnames(NB_Cell);
    for p = 1:numel(plateNames)
        meaPlate = plateNames{p};
        if ~isfield(NB_Cell.(meaPlate), DIV), continue; end
        wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));

        for w = 1:numel(wellNames)
            wellName = wellNames{w};
            NB = NB_Cell.(meaPlate).(DIV).(wellName);

            % --- Check well activity ---
            if ~isstruct(NB) || ~isfield(NB,'starts') || ...
                    isempty(NB.starts) || all(isnan(NB.starts))
                fprintf('[INACTIVE] %s | %s | %s | %s | %s — no NB detected, excluded.\n', ...
                    clName, detName, divName, meaPlate, wellName);
                continue
            end

            % --- Check numfragments field (correct field name) ---
            if ~isfield(NB,'num_fragments') || isempty(NB.num_fragments)
                fprintf('[NO_FRAGS_FIELD] %s | %s | %s | %s | %s — num_fragments missing, excluded.\n', ...
                    clName, detName, divName, meaPlate, wellName);
                continue
            end

            frags    = NB.num_fragments(:);
            frags    = frags(isfinite(frags));
            nTotal   = numel(frags);

            % --- Exclude NB with 0 fragments ---
            nSkipped = sum(frags == 0);
            if nSkipped > 0
                fprintf('[FRAG_SKIP] %s | %s | %s | %s | %s — %d/%d NBs have 0 fragments, excluded from ratio.\n', ...
                    clName, detName, divName, meaPlate, wellName, nSkipped, nTotal);
            end

            frags_valid = frags(frags >= 1);
            if isempty(frags_valid)
                fprintf('[NO_VALID_FRAGS] %s | %s | %s | %s | %s — no NB with frags >= 1, excluded.\n', ...
                    clName, detName, divName, meaPlate, wellName);
                continue
            end

            absVals(end+1) = sum(frags_valid);                        %#ok  Fig 1
            ratVals(end+1) = sum(frags_valid) / numel(frags_valid);   %#ok  Fig 2
        end
    end
end

function plotDivGrouped(DATA_T, DATA_D, CellLineNames, selectedDIVs, ...
        colorT, colorD, outlierK, yLabelStr, titleStr, figNameBase, subtitleStr, ...
        saveFig, saveFolder)

    nCL  = numel(CellLineNames);
    nDIV = numel(selectedDIVs);
    nCols = 2; nRows = ceil(nDIV / nCols);
    bw = 0.28; offT = -0.22; offD = 0.22;

    % --- Pairwise Tampere-vs-Doorn test per (cell line, DIV), FDR-corrected ---
    % (same statistical workflow as script3)
    pMat = nan(nCL, nDIV);
    for c = 1:nCL
        for d = 1:nDIV
            pMat(c,d) = run_pairwise_test(DATA_T{c,d}(:), DATA_D{c,d}(:), false);
        end
    end
    pAdjMat = reshape(fdr_bh(pMat(:)), nCL, nDIV);

    fig = figure('Color','w','Position',[60 60 1000 700]);
    tl  = tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');
%     title(tl, sprintf('%s (IQR k=%.1f)\n%s', titleStr, outlierK, subtitleStr), 'FontSize',11);

    for d = 1:nDIV
        ax = nexttile;
        hold(ax,'on');
        yAll = [];
        for c = 1:nCL
            dT = DATA_T{c,d}; dD = DATA_D{c,d};
            yAll = [yAll; dT(:); dD(:)]; %#ok
            if ~isempty(dT), drawBox(ax, c+offT, dT, colorT, bw); end
            if ~isempty(dD), drawBox(ax, c+offD, dD, colorD, bw); end
        end
        xticks(ax, 1:nCL); xticklabels(ax, CellLineNames);
        % Short cell-line names -> keep ticks horizontal so the slanted labels
        % of upper panels don't run into the panel below (compact tiling).
        ax.XAxis.FontSize = 10; xtickangle(ax, 0);
        xlim(ax, [0.5, nCL+0.5]);
        yspan = 0;
        if ~isempty(yAll) && any(isfinite(yAll))
            yspan = max(yAll(isfinite(yAll)));
            ylim(ax, [0, yspan*1.25]);
        end
        % --- Significance brackets (Tampere vs Doorn) above each cell line ---
        for c = 1:nCL
            pairVals = [DATA_T{c,d}(:); DATA_D{c,d}(:)];
            pairVals = pairVals(isfinite(pairVals));
            if isempty(pairVals) || isnan(pAdjMat(c,d)) || yspan==0, continue; end
            yBar = max(pairVals) + 0.06*yspan;
            drawSigStar(ax, c+offT, c+offD, yBar, 0.025*yspan, pAdjMat(c,d));
        end
%         title(ax, selectedDIVs{d}, 'FontSize',12,'FontWeight','bold');
        grid(ax,'on'); box(ax,'on');
        if d == 1
            hT = patch(ax,NaN,NaN,colorT,'EdgeColor',colorT*0.55,'LineWidth',1.2);
            hD = patch(ax,NaN,NaN,colorD,'EdgeColor',colorD*0.55,'LineWidth',1.2);
            legend(ax,[hT hD],{'Tampere','Doorn'},'Location','northeast','FontSize',10);
        end
    end
    for k = nDIV+1:nRows*nCols, ax = nexttile; axis(ax,'off'); end

    % Single shared y-axis label for the whole layout (avoids per-panel
    % ylabels colliding with the adjacent column in compact tiling).
    ylabel(tl, yLabelStr, 'FontSize',15,'FontWeight','bold');

    if saveFig
        exportgraphics(fig, fullfile(saveFolder,[figNameBase '.png']),'Resolution',300);
        fprintf('Saved: %s.png\n', figNameBase);
    end
end

function drawBox(ax, xCenter, data, clr, halfWidth)
data = data(isfinite(data));
if numel(data) < 2, return; end

q1 = prctile(data,25);
med = median(data);
q3 = prctile(data,75);
iqr_ = q3-q1;

wLow = max(data(data >= q1-1.5*iqr_));
wHig = min(data(data <= q3+1.5*iqr_));
outs = data(data < q1-1.5*iqr_ | data > q3+1.5*iqr_);

xl = xCenter-halfWidth;
xr = xCenter+halfWidth;
edgC = clr*0.55;

patch(ax,[xl xr xr xl xl],[q1 q1 q3 q3 q1],clr,...
    'EdgeColor',edgC,'LineWidth',1.3,'FaceAlpha',0.78);

line(ax,[xl xr],[med med],'Color',edgC,'LineWidth',2.0);
line(ax,[xCenter xCenter],[wLow q1],'Color',edgC,'LineWidth',1.2);
line(ax,[xCenter xCenter],[q3 wHig],'Color',edgC,'LineWidth',1.2);
line(ax,[xl xr],[wLow wLow],'Color',edgC,'LineWidth',1.0);
line(ax,[xl xr],[wHig wHig],'Color',edgC,'LineWidth',1.0);

m = mean(data);
plot(ax,xCenter,m,'d','Color',edgC,'MarkerFaceColor',edgC,...
    'MarkerSize',5,'HandleVisibility','off');

% CORRETTA:
dy = max([0.08*iqr_, 0.02*range(data), eps]);

if isempty(outs)
    topLimit = wHig;
else
    topLimit = max([wHig; outs(:)]);
end

if (topLimit - q3) >= dy
    yText = q3 + dy;
    vAlign = 'bottom';
else
    yText = q1 - dy;
    vAlign = 'top';
end

text(ax, xCenter, yText, sprintf('%.2f', m), ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment',vAlign, ...
    'FontSize',9, ...
    'FontWeight','bold', ...
    'Color',edgC, ...
    'BackgroundColor','w', ...
    'Margin',1, ...
    'Clipping','on');

if ~isempty(outs)
    plot(ax,xCenter*ones(size(outs)),outs,'.','Color',edgC,...
        'MarkerSize',8,'HandleVisibility','off');
end
end

function clean = removeOutliersIQR(data, k)
data = data(isfinite(data));
if numel(data) < 4
    clean = data;
    return
end
q1 = prctile(data, 25);
q3 = prctile(data, 75);
iqr_ = q3 - q1;
clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end

function drawSigStar(ax, xL, xR, yBar, tick, p_adj)
% Draw a significance bracket between xL and xR at height yBar, with the
% star string (or 'ns') just above it.
line(ax, [xL xL xR xR], [yBar-tick yBar yBar yBar-tick], ...
    'Color','k','LineWidth',1.1,'HandleVisibility','off');
s = starStr(p_adj);
if strcmp(s,'ns')
    text(ax, (xL+xR)/2, yBar+tick*0.5, s, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',9,'FontWeight','bold','Color',[0.45 0.45 0.45], ...
        'Clipping','on');
else
    text(ax, (xL+xR)/2, yBar+tick*0.2, s, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontSize',13,'FontWeight','bold','Color','k', ...
        'Clipping','on');
end
end

function s = starStr(p)
if isnan(p),        s = 'ns';
elseif p < 0.001,   s = '***';
elseif p < 0.01,    s = '**';
elseif p < 0.05,    s = '*';
else,               s = 'ns';
end
end

function [p, test_name, d_eff] = run_pairwise_test(xA, xB, paired)
% Adaptive two-group test (mirrors script3 workflow):
%   normality (Shapiro-Wilk -> Lilliefors fallback) decides
%   t-test / Welch (unpaired) or signed-rank / Mann-Whitney.
if nargin < 3, paired = false; end
test_name = ''; d_eff = NaN; p = NaN;
xA = xA(isfinite(xA)); xB = xB(isfinite(xB));
if numel(xA) < 2 || numel(xB) < 2, return; end

if paired
    n = min(numel(xA), numel(xB));
    xA = xA(1:n); xB = xB(1:n);
    diffs = xA - xB;
    if is_normal(diffs)
        [~, p] = ttest(xA, xB);  test_name = 'paired t-test';
    else
        p = signrank(xA, xB);    test_name = 'Wilcoxon signed-rank';
    end
    sd = std(diffs);
    if sd > 0, d_eff = mean(diffs)/sd; end
else
    normA = is_normal(xA); normB = is_normal(xB);
    if normA && normB
        [~, pv] = vartest2(xA, xB);
        if pv < 0.05
            [~, p] = ttest2(xA, xB, 'Vartype','unequal'); test_name = 'Welch t-test';
        else
            [~, p] = ttest2(xA, xB);                      test_name = 'Student t-test';
        end
    else
        p = ranksum(xA, xB); test_name = 'Mann-Whitney U';
    end
    nA = numel(xA); nB = numel(xB);
    sp = sqrt(((nA-1)*var(xA)+(nB-1)*var(xB))/max(nA+nB-2,1));
    if sp > 0, d_eff = (mean(xA)-mean(xB))/sp; end
end
end

function tf = is_normal(x)
x = x(isfinite(x));
if numel(x) < 3, tf = true; return; end
tf = true;
try
    [~, pN] = swtest(x);
    tf = pN > 0.05;
catch
    try
        [~, pN] = lillietest(x);
        tf = pN > 0.05;
    catch
        tf = true;
    end
end
end

function p_adj = fdr_bh(p)
% Benjamini-Hochberg FDR correction. NaNs are preserved.
p = p(:);
p_adj = nan(size(p));
valid = find(isfinite(p));
m = numel(valid);
if m == 0, return; end
[ps, order] = sort(p(valid));
ranks = (1:m)';
adj = ps .* m ./ ranks;
% enforce monotonicity (cumulative min from the largest)
for i = m-1:-1:1
    adj(i) = min(adj(i), adj(i+1));
end
adj = min(adj, 1);
out = nan(m,1);
out(order) = adj;
p_adj(valid) = out;
end
