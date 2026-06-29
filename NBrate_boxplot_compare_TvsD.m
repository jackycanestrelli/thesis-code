function NBrate_boxplot_compare_TvsD()
% NBrate_boxplot_compare_TvsD  NB-rate boxplots comparing Tampere vs Doorn detectors per DIV.
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

% BOXPLOT: NB Rate — Tampere vs Doorn
% Una figura per ogni DIV selezionato.
% Per ogni figura: cell lines sull'asse X, 2 box affiancati
%   Tampere (blue) | Doorn (red)
% 1 data point = 1 active well  |  NB rate = num_NB / 10 min
% Inactive wells (no NB) → excluded + [INACTIVE] debug message.
% REQUIRED IN WORKSPACE:
%   NB_Tampere_DD1C, NB_Tampere_DD5A, NB_Tampere_TUBA, NB_Tampere_C10902
%   NB_Doorn_DD1C,   NB_Doorn_DD5A,   NB_Doorn_TUBA,   NB_Doorn_C10902

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
recordingTime_s = 600;
recordingTime_m = recordingTime_s / 60;   % 10 min → NBs/min
outlierK        = 1.5;
saveFigures     = true;
saveFolder      = fullfile(outRoot, 'all_outputs', 'NBR_TvsD');

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
DATA_T_rate = cell(nCL, nDIV);
DATA_D_rate = cell(nCL, nDIV);

for c = 1:nCL
    n    = CellLineNames{c};
    NB_T = eval(sprintf('NB_Tampere_%s', n));
    NB_D = eval(sprintf('NB_Doorn_%s',   n));

    for d = 1:nDIV
        DIV = selectedDIVs{d};

        cntT = collectNBcount(NB_T, DIV, n, 'Tampere', DIV);
        cntD = collectNBcount(NB_D, DIV, n, 'Doorn',   DIV);

        if outlierK > 0
            cntT = removeOutliersIQR(cntT, outlierK);
            cntD = removeOutliersIQR(cntD, outlierK);
        end

        DATA_T_rate{c,d} = cntT / recordingTime_m;   % NBs/min
        DATA_D_rate{c,d} = cntD / recordingTime_m;

        fprintf('[DATA] %s | %s | Tampere: %d active wells | Doorn: %d active wells\n', ...
            n, DIV, numel(cntT), numel(cntD));
    end
end

%% ===== PAIRWISE STATS: Tampere vs Doorn (same workflow as script3) =====
% For each cell line and DIV, compare the per-well NB-rate distributions of the
% two detectors with the adaptive test used in script3 (normality -> t/Welch,
% otherwise Mann-Whitney), then BH-FDR correct across all comparisons of this
% metric. Stars are drawn above each Tampere/Doorn pair in the figures below.
pMat = nan(nCL, nDIV);
for c = 1:nCL
    for d = 1:nDIV
        pMat(c,d) = run_pairwise_test(DATA_T_rate{c,d}(:), DATA_D_rate{c,d}(:), false);
    end
end
pAdjMat = reshape(fdr_bh(pMat(:)), nCL, nDIV);

%% ===== ONE FIGURE PER DIV =====
bw   = 0.28;
offT = -0.22;
offD =  0.22;
yLabelStr = sprintf('NB rate (NBs/min)  [rec = %d s]', recordingTime_s);

for d = 1:nDIV
    DIV = selectedDIVs{d};

    fig = figure('Color','w','Position',[100 100 700 500]);
    ax  = axes(fig); %#ok
    hold(ax,'on');

    yAll = [];
    for c = 1:nCL
        dT = DATA_T_rate{c,d};
        dD = DATA_D_rate{c,d};
        yAll = [yAll; dT(:); dD(:)]; %#ok
        if ~isempty(dT), drawBox(ax, c+offT, dT, colorT, bw); end
        if ~isempty(dD), drawBox(ax, c+offD, dD, colorD, bw); end
    end

    xticks(ax, 1:nCL);
    xticklabels(ax, CellLineNames);
    xtickangle(ax, 0);   % short names -> horizontal, avoids x-axis label overlap
    ax.XAxis.FontSize = 11;
    xlim(ax, [0.5, nCL+0.5]);
    yspan = 0;
    if ~isempty(yAll) && any(isfinite(yAll))
        yspan = max(yAll(isfinite(yAll)));
        ylim(ax, [0, yspan*1.25]);   % extra headroom for significance stars
    end

    % --- significance brackets: Tampere vs Doorn, one per cell line ---
    for c = 1:nCL
        pairVals = [DATA_T_rate{c,d}(:); DATA_D_rate{c,d}(:)];
        pairVals = pairVals(isfinite(pairVals));
        if isempty(pairVals) || isnan(pAdjMat(c,d)) || yspan==0, continue; end
        yBar = max(pairVals) + 0.06*yspan;
        drawSigStar(ax, c+offT, c+offD, yBar, 0.025*yspan, pAdjMat(c,d));
    end

    ylabel(ax, yLabelStr, 'FontSize',15);
%     title(ax, sprintf('NB Rate — Tampere vs Doorn  |  %s', DIV), ...
%         'FontSize',13,'FontWeight','bold');
    subtitle(ax, sprintf('IQR k=%.1f  |  1 data point = 1 well  |  num\\_NB / 10 min', outlierK), ...
        'FontSize',9,'Color',[0.4 0.4 0.4]);
    grid(ax,'on'); box(ax,'on');

    hT = patch(ax,NaN,NaN,colorT,'EdgeColor',colorT*0.55,'LineWidth',1.2);
    hD = patch(ax,NaN,NaN,colorD,'EdgeColor',colorD*0.55,'LineWidth',1.2);
    legend(ax,[hT hD],{'Tampere','Doorn'},'Location','northeast','FontSize',11);

    if saveFigures
        fname = sprintf('boxplot_NBrate_TampereVsDoorn_%s.png', DIV);
        exportgraphics(fig, fullfile(saveFolder, fname), 'Resolution',300);
        fprintf('Saved: %s\n', fname);
    end
end

% LOCAL FUNCTIONS

end

function vals = collectNBcount(NB_Cell, DIV, clName, detName, divName)
    vals = [];
    if isempty(NB_Cell), return; end
    plateNames = fieldnames(NB_Cell);
    for p = 1:numel(plateNames)
        meaPlate = plateNames{p};
        if ~isfield(NB_Cell.(meaPlate), DIV), continue; end
        wellNames = fieldnames(NB_Cell.(meaPlate).(DIV));
        for w = 1:numel(wellNames)
            wellName = wellNames{w};
            NB = NB_Cell.(meaPlate).(DIV).(wellName);

            nNB = 0;
            if isstruct(NB) && isfield(NB,'starts') && ~isempty(NB.starts) && ~all(isnan(NB.starts))
                nNB = numel(NB.starts);
            elseif isstruct(NB) && isfield(NB,'num_bursts') && ~isempty(NB.num_bursts)
                nNB = NB.num_bursts;
            elseif isnumeric(NB) && size(NB,2) >= 2
                nNB = size(NB,1);
            end

            if nNB == 0
                fprintf('[INACTIVE] %s | %s | %s | %s | %s — no NB detected, excluded.\n', ...
                    clName, detName, divName, meaPlate, wellName);
                continue
            end
            vals(end+1) = nNB; %#ok
        end
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
    if numel(data) < 4, clean = data; return; end
    q1 = prctile(data,25); q3 = prctile(data,75); iqr_ = q3-q1;
    clean = data(data >= q1-k*iqr_ & data <= q3+k*iqr_);
end

% ===== Stats helpers (mirror script3) =====================================
function drawSigStar(ax, xL, xR, yBar, tick, p_adj)
% Bracket from xL to xR at height yBar with the significance star above it.
s = starStr(p_adj);
if s == "ns"
    barCol = [0.55 0.55 0.55]; fs = 9;
else
    barCol = [0.10 0.10 0.10]; fs = 13;
end
line(ax, [xL xL], [yBar-tick yBar], 'Color',barCol,'LineWidth',1.1,'HandleVisibility','off');
line(ax, [xR xR], [yBar-tick yBar], 'Color',barCol,'LineWidth',1.1,'HandleVisibility','off');
line(ax, [xL xR], [yBar yBar],      'Color',barCol,'LineWidth',1.1,'HandleVisibility','off');
text(ax, (xL+xR)/2, yBar+0.2*tick, char(s), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'FontSize',fs,'FontWeight','bold','Color',barCol,'Clipping','on');
end

function s = starStr(p)
if     isnan(p);  s = "";
elseif p<0.001;   s = "***";
elseif p<0.01;    s = "**";
elseif p<0.05;    s = "*";
else;             s = "ns";
end
end

function [p, test_name, d_eff] = run_pairwise_test(xA, xB, paired)
xA = xA(~isnan(xA)); xB = xB(~isnan(xB));
p = NaN; test_name = 'insufficient_n'; d_eff = NaN;
if numel(xA) < 3 || numel(xB) < 3; return; end
nA_ok = is_normal(xA); nB_ok = is_normal(xB);
if nA_ok && nB_ok
    try; equal_var = vartest2(xA,xB) > 0.05; catch; equal_var = false; end
    if paired;        [~,p] = ttest(xA,xB);                       test_name = 'paired t-test';
    elseif equal_var; [~,p] = ttest2(xA,xB);                      test_name = 't-test (equal var)';
    else;             [~,p] = ttest2(xA,xB,'Vartype','unequal');  test_name = 'Welch t-test'; end
else
    if paired
        try; p = signrank(xA,xB); test_name = 'Wilcoxon signed-rank';
        catch; p = ranksum(xA,xB); test_name = 'Mann-Whitney (fallback)'; end
    else; p = ranksum(xA,xB); test_name = 'Mann-Whitney U'; end
end
sA = std(xA,'omitnan'); sB = std(xB,'omitnan');
nA = numel(xA); nB = numel(xB);
sp = sqrt(((nA-1)*sA^2+(nB-1)*sB^2)/max(nA+nB-2,1));
if sp > 0; d_eff = (mean(xA,'omitnan')-mean(xB,'omitnan'))/sp; end
end

function tf = is_normal(x)
x = x(~isnan(x));
if numel(x) < 4; tf = false; return; end
try; [~,p] = swtest(x,0.05);
catch; try; [~,p] = lillietest(x); catch; tf=false; return; end; end
tf = p > 0.05;
end

function p_adj = fdr_bh(p)
p = p(:); n = numel(p);
valid = ~isnan(p); pv = p(valid);
[ps,si] = sort(pv); m = numel(ps);
if m==0; p_adj=p; return; end
p_adj_v = nan(m,1);
for i = m:-1:1
    if i==m; p_adj_v(i) = ps(i);
    else;    p_adj_v(i) = min(p_adj_v(i+1), ps(i)*m/i); end
end
p_adj_v = min(p_adj_v,1);
p_adj_sorted = nan(m,1); p_adj_sorted(si) = p_adj_v;
p_adj = nan(n,1); p_adj(valid) = p_adj_sorted;
end
