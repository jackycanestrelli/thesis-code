function stats_RNDSPIKES_new()
% stats_RNDSPIKES_new  Pairwise statistics (LME/ANOVA fallback) and heatmaps for % random spikes.
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
reqVars = {'DD1C','DD5A','TUBA','C10902','SampFreq','ParamsBurst','recDuration'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% STAT HEATMAP: % Random Spikes — 4 Cell Lines pairwise
% Approccio: identico a stats_analysis_def.m (Analisi A1)
% - LME overall (CL + DIV effect), fallback 2-way ANOVA
% - Post-hoc pairwise con test adattivo (normalità → t/Welch, altrimenti Mann-Whitney)
% - FDR (Benjamini-Hochberg) correction
% - Output: heatmap -log10(p_adj) | bianco=ns, grigio=NA, blu scuro=sig
% REQUIRED IN WORKSPACE:
% DD1C, DD5A, TUBA, C10902
% ParamsBurst, SampFreq, recDuration

%% ===== USER SETTINGS =====
selectedDIVs  = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
DIVsnumplot   = [53 60 67 83 90];
outlierK      = 1.5;
recordingTime = 600;
saveFigure    = true;
saveFolder    = fullfile(outRoot, 'all_outputs', 'rndspikes_stats');

if ~exist('recDuration','var') || isempty(recDuration)
    warning('recDuration non trovato. Uso recordingTime = %d s.', recordingTime);
    recDuration = recordingTime;
end

if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

cellLineNames = {'DD1C','DD5A','TUBA','C10902'};
nCL  = numel(cellLineNames);
nDIV = numel(selectedDIVs);

%% ===== BUILD CELL LINE LIST =====
CellLines = {};
for i = 1:nCL
    n = cellLineNames{i};
    if exist(n,'var')
        CellLines{end+1} = struct('name', n, 'Data', eval(n)); %#ok
    else
        fprintf('[SKIP] %s not found in workspace\n', n);
    end
end
nCLvalid = numel(CellLines);
if nCLvalid == 0, error('No cell lines found'); end

%% ===== DATA COLLECTION =====
dataM = cell(nCLvalid, nDIV);

for c = 1:nCLvalid
    cellData = CellLines{c}.Data;
    clName   = CellLines{c}.name;
    plateNames = fieldnames(cellData);

    for d = 1:nDIV
        DIV  = selectedDIVs{d};
        vals = [];

        for p = 1:numel(plateNames)
            pl = plateNames{p};
            if ~isfield(cellData.(pl), DIV), continue; end
            wellNames = fieldnames(cellData.(pl).(DIV));

            for w = 1:numel(wellNames)
                spikeTable = cellData.(pl).(DIV).(wellNames{w});
                if ~istable(spikeTable) || height(spikeTable)==0 || ...
                   ~ismember('Channel', spikeTable.Properties.VariableNames) || ...
                   ~ismember('Time',    spikeTable.Properties.VariableNames)
                    continue
                end
                Time    = spikeTable.Time;
                Channel = spikeTable.Channel;
                valid   = isfinite(Time) & Time >= 0;
                Time    = Time(valid); Channel = Channel(valid);
                if isempty(Time), continue; end

                pct = getMeanRandomSpikePercentageWell(Time, Channel, SampFreq, ParamsBurst, recDuration);
                if isnan(pct), continue; end
                vals(end+1) = pct; %#ok
            end
        end

        if outlierK > 0 && numel(vals) >= 4
            vals = removeOutliersIQR(vals, outlierK);
        end
        dataM{c,d} = vals;
        fprintf('[DATA] %s | %s | n=%d\n', clName, DIV, numel(vals));
    end
end

%% ===== OVERALL MODEL (LME / ANOVA) =====
allVals = []; clVec = []; divVec = []; subjVec = [];
subjID  = 0;
for c = 1:nCLvalid
    for d = 1:nDIV
        v = dataM{c,d};
        if isempty(v), continue; end
        n = numel(v);
        allVals  = [allVals;  v(:)];               %#ok
        clVec    = [clVec;    repmat(c,n,1)];       %#ok
        divVec   = [divVec;   repmat(d,n,1)];       %#ok
        subjVec  = [subjVec;  subjID+(1:n)'];       %#ok
        subjID   = subjID + n;
    end
end

valid   = ~isnan(allVals);
allVals = allVals(valid); clVec = clVec(valid);
divVec  = divVec(valid);  subjVec = subjVec(valid);

CLcat   = categorical(clVec);
DIVcat  = categorical(divVec);
Subjcat = categorical(subjVec);
tbl     = table(allVals, CLcat, DIVcat, Subjcat, ...
                'VariableNames',{'Val','CL','DIV','Subj'});

try
    lme       = fitlme(tbl, 'Val ~ CL + DIV + (1|Subj)', 'FitMethod','REML');
    at        = anova(lme);
    pcl       = at.pValue(2);
    modelName = 'Mixed effects';
catch ME
    warning('fitlme failed: %s — fallback ANOVA', ME.message);
    [pan,~,~] = anovan(allVals, {clVec, divVec}, 'model','interaction', ...
                       'varnames',{'CL','DIV'}, 'display','off');
    pcl       = pan(1);
    modelName = '2-way ANOVA fallback';
end
fprintf('\nOverall CL effect: p=%.4f %s [%s]\n', pcl, getStar(pcl), modelName);

%% ===== PAIRWISE TEST + FDR =====
pairIdx    = nchoosek(1:nCLvalid, 2);
nPairs     = size(pairIdx, 1);

pairLabels = cell(1, nPairs);
for pp = 1:nPairs
    pairLabels{pp} = sprintf('%s vs %s', ...
        CellLines{pairIdx(pp,1)}.name, CellLines{pairIdx(pp,2)}.name);
end

% Collect raw p-values and metadata
praw    = nan(nDIV, nPairs);
signmat = zeros(nDIV, nPairs);
hasData = false(nDIV, nPairs);
testNames = cell(nDIV, nPairs);
cohenD  = nan(nDIV, nPairs);
meanDiff= nan(nDIV, nPairs);

for d = 1:nDIV
    for pp = 1:nPairs
        g1 = dataM{pairIdx(pp,1), d};
        g2 = dataM{pairIdx(pp,2), d};
        if numel(g1) >= 3 && numel(g2) >= 3
            hasData(d,pp)  = true;
            signmat(d,pp)  = sign(median(g1) - median(g2));
            [p_, tname_, d_] = run_pairwise_test(g1, g2, false);
            praw(d,pp)     = p_;
            testNames{d,pp}= tname_;
            cohenD(d,pp)   = d_;
            meanDiff(d,pp) = mean(g1,'omitnan') - mean(g2,'omitnan');
        end
    end
end

% FDR across all valid pairs (all DIVs together)
praw_vec = praw(:);
padj_vec = fdr_bh(praw_vec);
padj     = reshape(padj_vec, nDIV, nPairs);

%% ===== SAVE PAIRWISE TABLE =====
rows = {};
for d = 1:nDIV
    for pp = 1:nPairs
        if hasData(d,pp)
            rows{end+1} = {selectedDIVs{d}, pairLabels{pp}, ...
                praw(d,pp), padj(d,pp), getStar(padj(d,pp)), ...
                cohenD(d,pp), meanDiff(d,pp), testNames{d,pp}}; %#ok
        end
    end
end
if ~isempty(rows)
    Tpair = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'DIV','Pair','p_raw','p_adj_FDR','Significance','Cohen_d','MeanDiff','Test'});
    writetable(Tpair, fullfile(saveFolder,'RndSpikes_pairwise.xlsx'));
    fprintf('Saved: RndSpikes_pairwise.xlsx\n');
end

%% ===== PLOT HEATMAP =====
% Build matrix M = -log10(p_adj) for significant cells, NaN for no data
M      = nan(nDIV, nPairs);
Lstars = strings(nDIV, nPairs);
Ldir   = strings(nDIV, nPairs);

for d = 1:nDIV
    for pp = 1:nPairs
        if ~hasData(d,pp), continue; end  % stays NaN = grigio
        p_ = padj(d,pp);
        if isnan(p_) || p_ >= 0.05
            M(d,pp) = 0;                  % ns = bianco (log10(1)=0)
        else
            M(d,pp) = min(-log10(p_), 6);
            cl1 = pairIdx(pp,1); cl2 = pairIdx(pp,2);
            Lstars(d,pp) = getStar(p_);
            if signmat(d,pp) > 0
                Ldir(d,pp) = sprintf('%s>%s', CellLines{cl1}.name, CellLines{cl2}.name);
            else
                Ldir(d,pp) = sprintf('%s>%s', CellLines{cl2}.name, CellLines{cl1}.name);
            end
        end
    end
end

plot_logp_heatmap(M, Lstars, Ldir, ...
    pairLabels, ...
    arrayfun(@(x) sprintf('DIV%d',x), DIVsnumplot, 'UniformOutput',false), ...
    sprintf('%%Random Spikes — Pairwise comparison (FDR BH)\n[%s]', modelName), ...
    'Cell line pair', 'DIV', ...
    fullfile(saveFolder,'heatmap_pairwise_RndSpikes'), ...
    saveFigure);

%% LOCAL FUNCTIONS

end

function [p, test_name, d_eff] = run_pairwise_test(xA, xB, paired)
xA = xA(~isnan(xA)); xB = xB(~isnan(xB));
p = NaN; test_name = 'insufficient_n'; d_eff = NaN;
if numel(xA) < 3 || numel(xB) < 3; return; end

nA_ok = is_normal(xA);
nB_ok = is_normal(xB);

if nA_ok && nB_ok
    try; equal_var = vartest2(xA, xB) > 0.05; catch; equal_var = false; end
    if paired
        [~, p] = ttest(xA, xB); test_name = 'paired t-test';
    elseif equal_var
        [~, p] = ttest2(xA, xB); test_name = 't-test (equal var)';
    else
        [~, p] = ttest2(xA, xB, 'Vartype','unequal'); test_name = 'Welch t-test';
    end
else
    if paired
        try; p = signrank(xA,xB); test_name = 'Wilcoxon signed-rank';
        catch; p = ranksum(xA,xB); test_name = 'Mann-Whitney (fallback)'; end
    else
        p = ranksum(xA, xB); test_name = 'Mann-Whitney U';
    end
end

sA = std(xA,'omitnan'); sB = std(xB,'omitnan');
nA = numel(xA); nB = numel(xB);
sp = sqrt(((nA-1)*sA^2 + (nB-1)*sB^2) / max(nA+nB-2,1));
if sp > 0; d_eff = (mean(xA,'omitnan') - mean(xB,'omitnan')) / sp; end
end

function tf = is_normal(x)
x = x(~isnan(x));
if numel(x) < 4; tf = false; return; end
try; [~,p] = swtest(x, 0.05);
catch; try; [~,p] = lillietest(x); catch; tf = false; return; end; end
tf = p > 0.05;
end

function p_adj = fdr_bh(p)
p = p(:); n = numel(p);
valid = ~isnan(p); pv = p(valid);
[ps, si] = sort(pv); m = numel(ps);
if m == 0; p_adj = p; return; end
p_adj_v = nan(m,1);
for i = m:-1:1
    if i == m; p_adj_v(i) = ps(i);
    else;      p_adj_v(i) = min(p_adj_v(i+1), ps(i)*m/i); end
end
p_adj_v = min(p_adj_v,1);
p_adj_sorted = nan(m,1); p_adj_sorted(si) = p_adj_v;
p_adj = nan(n,1); p_adj(valid) = p_adj_sorted;
end

function s = getStar(p)
s = '';
if isnan(p) || isempty(p), return; end
if     p < 0.001, s = '****';  % kept for legacy; BH rarely gives 0.0001
elseif p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
end
end

function plot_logp_heatmap(M, Lstars, Ldir, xLabels, yLabels, ttl, xlbl, ylbl, savebase, doSave)
% M(row,col): NaN=no data(grigio), 0=ns(bianco), >0=-log10(p_adj) scuro
[nrows, ncols] = size(M);
nColors = 256;

% Build colormap: 1=grigio(NA), 2=bianco(ns), 3..256=blu crescente
cmap = zeros(nColors, 3);
cmap(1,:) = [0.75 0.75 0.75];          % grigio = no data
cmap(2,:) = [1.00 1.00 1.00];          % bianco = ns
blueStart  = [0.90 0.93 1.00];
blueEnd    = [0.05 0.18 0.65];
for i = 1:nColors-2
    frac = (i-1)/(nColors-3);
    cmap(i+2,:) = blueStart + frac*(blueEnd - blueStart);
end

maxVal  = 6;
imgData = ones(nrows, ncols);          % default = 1 (grigio)

for r = 1:nrows
    for c = 1:ncols
        if isnan(M(r,c))
            imgData(r,c) = 1;          % grigio = no data
        elseif M(r,c) == 0
            imgData(r,c) = 2;          % bianco = ns
        else
            frac = min(M(r,c), maxVal) / maxVal;
            imgData(r,c) = round(3 + frac*(nColors-3));
            imgData(r,c) = min(max(imgData(r,c), 3), nColors);
        end
    end
end

figW = max(400, 130*ncols + 120);
figH = max(300, 90*nrows  + 120);
fig  = figure('Color','w','Position',[100 100 figW figH]);
ax   = axes('Parent', fig);
imagesc(ax, imgData);
colormap(ax, cmap);
clim(ax, [1 nColors]);
hold(ax,'on');

% Grid lines
for r = 0.5:1:nrows+0.5
    line(ax,[0.5 ncols+0.5],[r r],'Color','k','LineWidth',0.7,'HandleVisibility','off');
end
for c = 0.5:1:ncols+0.5
    line(ax,[c c],[0.5 nrows+0.5],'Color','k','LineWidth',0.7,'HandleVisibility','off');
end

% NOTE: no text is drawn inside the cells (no stars / direction / nd / ns).
% Significance is encoded by the colour scale (-log10 p_adj); "ns" cells stay
% white and "nd" (not-determined) cells stay grey, as defined by the colormap.
% This keeps the heatmap clean and avoids text overlapping in small cells.

% Colorbar
cb = colorbar(ax);
nTicks = 7;
cb.Ticks      = linspace(2, nColors, nTicks);
cb.TickLabels = {'0 (ns)', '1', '2', '3', '4', '5', sprintf('>=%d', maxVal)};
cb.Label.String   = '-log_{10}(p_{adj})';
cb.Label.FontSize = 11;

set(ax, 'XTick', 1:ncols, 'XTickLabel', xLabels, 'XTickLabelRotation', 35, ...
        'YTick', 1:nrows,  'YTickLabel', yLabels, ...
        'FontSize', 10, 'TickLength', [0 0], 'LineWidth', 0.8);
xlabel(ax, xlbl, 'FontSize',15, 'FontWeight','bold');
ylabel(ax, ylbl, 'FontSize',15, 'FontWeight','bold');
% title(ax, ttl,   'FontSize',12, 'FontWeight','bold');
box(ax,'on'); set(ax,'Layer','top');

drawnow;
if doSave
    exportgraphics(fig, [savebase '.png'], 'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'], 'ContentType','vector');
    fprintf('Saved: %s\n', savebase);
end
end

function pctRandom = getMeanRandomSpikePercentageWell(Time, Channel, SampFreq, ParamsBurst, recDuration)
Channel = string(Channel);
electrodes = unique(Channel);
totalSpikesAll = 0; burstSpikesAll = 0;
t0 = min(Time);
for e = 1:numel(electrodes)
    elec = electrodes(e);
    spT  = sort(Time(Channel == elec) - t0);
    sps  = unique(round(spT * SampFreq) + 1);
    sps  = sps(sps >= 1);
    if isempty(sps), continue; end
    totalSpikesAll = totalSpikesAll + numel(sps);
    if numel(sps) < 3, continue; end
    bt = DetectBurstsNeuro(sps, SampFreq, recDuration, ParamsBurst);
    if isempty(bt) || size(bt,1)==0, continue; end
    inBurst = false(size(sps));
    for b = 1:size(bt,1)
        inBurst = inBurst | (sps >= bt(b,1) & sps <= bt(b,2));
    end
    burstSpikesAll = burstSpikesAll + sum(inBurst);
end
if totalSpikesAll == 0; pctRandom = NaN;
else; pctRandom = ((totalSpikesAll - burstSpikesAll) / totalSpikesAll) * 100;
end
end

function clean = removeOutliersIQR(data, k)
data = data(~isnan(data));
if numel(data) < 4; clean = data; return; end
q1 = prctile(data,25); q3 = prctile(data,75); iqr_ = q3-q1;
clean = data(data >= q1-k*iqr_ & data <= q3+k*iqr_);
end
