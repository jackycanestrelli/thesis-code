function outputsection4_LME2()
% outputsection4_LME2  Per-cell-line boxplots comparing the four parameter sets
%   (X = metrics, colour = Params1..4) with an LME statistical workflow aligned
%   to script3: omnibus LME (metric ~ ParamSet + (1|Plate); Kruskal-Wallis
%   fallback when n<6), adaptive pairwise post-hoc (t/Welch or Mann-Whitney) with
%   Benjamini-Hochberg FDR and Cohen's d. Statistics use the FULL dataset; the
%   boxplots use IQR-trimmed data (k=1.5) for display only.
%
% INPUT : reads <outRoot>/Params1..4/metrics_per_well.csv (must include a Plate
%         column). Output root set via RESULTS_ROOT in the base workspace
%         (default <pwd>/results).
% OUTPUT: one boxplot per cell line (PNG + FIG) plus statistical_results_summary.csv
%         and posthoc_results.csv, under <outRoot>/Section 4 def/Boxplots_ByCellLine.

if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

Folder       = fullfile(outRoot, 'Section 4 def');
% Root holding the per-parameter-set metric CSVs (paramsX/metrics_per_well.csv).
% Override by setting METRICS_ROOT in the base workspace if the drive changes.
if evalin('base', "exist('METRICS_ROOT','var')")
    baseFolder = evalin('base','METRICS_ROOT');
else
    baseFolder = 'S:\met_narkilahti_neuro_sto-3700\MEA_data_internship\Jacopo\TEST SCRIPT APP';
end
paramFolders = {'params1','params2','params3','params4'};

selectedDIVs  = [53, 62, 67, 83, 88];
cellLineNames = {'DD1C','DD5A','TUBA','C10902'};
outlierK      = 1.5;     % box display only, not the statistics

outFolder = fullfile(outRoot, 'all_outputs', 'section4_boxplots');
if ~exist(outFolder,'dir'), mkdir(outFolder); end

plotJobs = {
    'NBR',       'TampereV1', 'NBR Tampere';
    'NBD',       'TampereV1', 'NBD';
    'CV_IBI',    'TampereV1', 'CV IBI';
    'FragPerNB', 'TampereV1', 'Frag per NB';
};

nMetrics = size(plotJobs,1);
nParams  = numel(paramFolders);

paramColors = [0.31 0.58 0.80;
               0.93 0.55 0.29;
               0.42 0.68 0.31;
               0.67 0.44 0.70];

offsets = linspace(-0.27, 0.27, nParams);

allData = table();
for i = 1:nParams
    csvPath = fullfile(baseFolder, paramFolders{i}, 'metrics_per_well.csv');
    if ~isfile(csvPath)
        warning('File non trovato: %s', csvPath); continue;
    end
    Ttemp = readtable(csvPath, 'TextType','string');
    Ttemp.ParamSet = repmat(string(paramFolders{i}), height(Ttemp), 1);
    allData = [allData; Ttemp]; %#ok
end

if isempty(allData), error('Nessun dato caricato.'); end

allData.CellLine = string(allData.CellLine);
allData.Detector = string(allData.Detector);
allData.ParamSet = string(allData.ParamSet);
allData.DIV      = double(allData.DIV);
if ~ismember('Plate', allData.Properties.VariableNames)
    error(['Colonna "Plate" assente nei CSV: serve per il random intercept (1|Plate). ' ...
           'Rigenera i metrics_per_well.csv con la colonna Plate, come in script3.']);
end
allData.Plate = string(allData.Plate);
allData = allData(ismember(allData.DIV, selectedDIVs), :);

resultsAll     = table();
resultsPosthoc = table();

for cl = 1:numel(cellLineNames)
    CL  = cellLineNames{cl};
    Tcl = allData(allData.CellLine == CL, :);

    dataMatrix = cell(nMetrics, nParams);
    fullMatrix = cell(nMetrics, nParams);
    fullPlate  = cell(nMetrics, nParams);
    for m = 1:nMetrics
        metricCol    = plotJobs{m,1};
        detectorName = plotJobs{m,2};
        Tsub = Tcl(Tcl.Detector == detectorName, :);
        Y = Tsub.(metricCol);
        if isstring(Y)||iscell(Y), Y = str2double(string(Y)); end
        Tsub.(metricCol) = Y;
        for p = 1:nParams
            idx  = Tsub.ParamSet == paramFolders{p};
            vals = Y(idx);
            pl   = Tsub.Plate(idx);
            ok   = isfinite(vals);
            fullMatrix{m,p} = vals(ok);                       % full data (statistics)
            fullPlate{m,p}  = pl(ok);
            dataMatrix{m,p} = removeOutliersIQR(vals(ok), outlierK);  % IQR-trimmed (display)
        end
    end

    allValsFlat = [];
    for m = 1:nMetrics
        for p = 1:nParams
            allValsFlat = [allValsFlat; dataMatrix{m,p}]; %#ok
        end
    end
    yMaxCL = max(allValsFlat,[],'omitnan');
    yMinCL = min(allValsFlat,[],'omitnan');
    if isnan(yMaxCL), yMaxCL = 1; end
    if isnan(yMinCL), yMinCL = 0; end
    yRange = yMaxCL - yMinCL;
    if yRange == 0, yRange = max(1, abs(yMaxCL)*0.2+1); end

    fig = figure('Color','w','Position',[100 100 1200 580]);
    ax  = axes(fig);
    hold(ax,'on');

    legendHandles = gobjects(0);
    legendNames   = strings(0);
    alreadyAdded  = false(1, nParams);
    wHighMat      = NaN(nMetrics, nParams);

    for m = 1:nMetrics
        xCenter = m;
        for p = 1:nParams
            vals = dataMatrix{m,p};
            if isempty(vals), continue; end

            xPos = xCenter + offsets(p);
            clr  = paramColors(p,:);

            drawStyledBox(ax, xPos, vals, clr, 0.13);

            q3v  = prctile(vals,75);
            iqrV = iqr(vals);
            wH   = min(max(vals(vals <= q3v + 1.5*iqrV)), q3v + 1.5*iqrV);
            wHighMat(m,p) = wH;

            text(ax, xPos, wH + 0.03*yRange, ...
                sprintf('%.2f', mean(vals,'omitnan')), ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment',  'bottom', ...
                'FontSize', 7.5, 'Color', clr*0.70);

            if ~alreadyAdded(p)
                h = patch(ax, NaN, NaN, clr, ...
                    'FaceAlpha',0.45,'EdgeColor',clr,'LineWidth',1.2);
                legendHandles(end+1) = h;                        %#ok
                legendNames(end+1)   = string(paramFolders{p}); %#ok
                alreadyAdded(p)      = true;
            end
        end
    end

    yBarBase = yMaxCL + 0.15*yRange;
    yBarStep = 0.07*yRange;
    pairs    = nchoosek(1:nParams, 2);
    nPairs   = size(pairs,1);

    for m = 1:nMetrics
        xCenter   = m;
        metricTag = sprintf('%s_%s', plotJobs{m,1}, plotJobs{m,2});

        yFull = []; gFull = strings(0); pFull = strings(0);
        for p = 1:nParams
            v = fullMatrix{m,p};
            if ~isempty(v)
                yFull = [yFull; v];                                          %#ok
                gFull = [gFull; repmat(string(paramFolders{p}),numel(v),1)]; %#ok
                pFull = [pFull; string(fullPlate{m,p})];                     %#ok
            end
        end
        if numel(unique(gFull)) < 2, continue; end

        Tlme = table(yFull, gFull, pFull, ...
            'VariableNames',{'y','ParamSet','Plate'});
        [pGlobal, test_used] = run_overall_lme_factor(Tlme, 'y', 'ParamSet', 'Plate', paramFolders);

        fprintf('\n%s | %s | Global: %s  p=%.4g\n', CL, metricTag, test_used, pGlobal);

        starGlobal = pValueToStar(pGlobal);
        wHighMax   = max(wHighMat(m,:),[],'omitnan');
        text(ax, xCenter, wHighMax + 0.10*yRange, starGlobal, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment',  'bottom', ...
            'FontSize', 13, 'FontWeight','bold', 'Color','k');

        resultsAll = [resultsAll; table( ...
            string(CL), string(metricTag), pGlobal, string(test_used), ...
            'VariableNames',{'CellLine','Metric','p_global','Test_global'})]; %#ok

        if pGlobal < 0.05
            pRaw      = nan(nPairs,1);
            testNames = strings(nPairs,1);
            dEff      = nan(nPairs,1);

            for rr = 1:nPairs
                vA = fullMatrix{m, pairs(rr,1)};
                vB = fullMatrix{m, pairs(rr,2)};
                [pRaw(rr), testNames(rr), dEff(rr)] = run_pairwise_test(vA, vB, false);
            end

            pAdj = fdr_bh(pRaw);

            barLevel = 0;
            for rr = 1:nPairs
                p1       = pairs(rr,1);
                p2       = pairs(rr,2);
                starPost = pValueToStar(pAdj(rr));

                fprintf('  %s vs %s | %s | p_raw=%.4g  p_adj=%.4g  %s\n', ...
                    paramFolders{p1}, paramFolders{p2}, ...
                    testNames(rr), pRaw(rr), pAdj(rr), starPost);

                resultsPosthoc = [resultsPosthoc; table( ...
                    string(CL), string(metricTag), ...
                    string(paramFolders{p1}), string(paramFolders{p2}), ...
                    testNames(rr), pRaw(rr), pAdj(rr), dEff(rr), string(starPost), ...
                    'VariableNames', {'CellLine','Metric','ParamA','ParamB', ...
                                      'Test','p_raw','p_adj_FDR','Cohen_d','Significance'})]; %#ok

                x1 = xCenter + offsets(p1);
                x2 = xCenter + offsets(p2);
                yB = yBarBase + barLevel * yBarStep;
                barLevel = barLevel + 1;
                capH = 0.03*yRange;

                cA   = paramColors(p1,:);
                cB   = paramColors(p2,:);
                isns = strcmp(starPost,'ns');
                if isns; barCol = [0.70 0.70 0.70]; legW = 1.0;
                else;    barCol = [0.15 0.15 0.15]; legW = 1.8; end
                line(ax, [x1,x2], [yB,yB],      'Color',barCol,'LineWidth',1.0); % connector
                line(ax, [x1,x1], [yB-capH,yB], 'Color',cA,'LineWidth',legW);    % gamba set A
                line(ax, [x2,x2], [yB-capH,yB], 'Color',cB,'LineWidth',legW);    % gamba set B
                if isns
                    text(ax, (x1+x2)/2, yB+0.005*yRange, 'ns', ...
                        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                        'FontSize', 8, 'Color', [0.55 0.55 0.55]);
                else
                    text(ax, (x1+x2)/2, yB+0.005*yRange, starPost, ...
                        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                        'FontSize', 11, 'FontWeight','bold', 'Color','k');
                end
            end
        end
    end

    nMaxBars   = nPairs;
    yTopNeeded = yBarBase + nMaxBars*yBarStep + 0.05*yRange;

    ylim(ax, [max(0, yMinCL - 0.05*yRange), yTopNeeded]);
    xlim(ax, [0.5, nMetrics + 0.5]);
    xticks(ax, 1:nMetrics);
    xticklabels(ax, plotJobs(:,3));
    xtickangle(ax, 0);
    grid(ax,'on');
    box(ax,'on');
    set(ax,'FontSize',14,'LineWidth',1,'GridAlpha',0.25);

    xlabel(ax, 'Metric', 'FontWeight','bold');
    ylabel(ax, 'Value',  'FontWeight','bold');
%     title(ax, sprintf('%s — comparison of 4 sets of parameters  (IQR k=%.1f)', CL, outlierK), ...
%         'FontWeight','bold','Interpreter','none');

    if ~isempty(legendHandles)
        legend(ax, legendHandles, cellstr(legendNames), ...
            'Location','northeastoutside','FontSize',9);
    end
    hold(ax,'off');

    fname = sprintf('%s_metrics_boxplot.png', CL);
    saveas(fig,  fullfile(outFolder, fname));
    savefig(fig, fullfile(outFolder, strrep(fname,'.png','.fig')));
    fprintf('Salvata: %s\n', fname);

end % cell-line loop

if ~isempty(resultsAll)
    writetable(resultsAll, fullfile(outFolder,'statistical_results_summary.csv'));
    fprintf('\nSalvato: statistical_results_summary.csv  (%d righe)\n', height(resultsAll));
end

if ~isempty(resultsPosthoc)
    writetable(resultsPosthoc, fullfile(outFolder,'posthoc_results.csv'));
    fprintf('Salvato: posthoc_results.csv  (%d righe)\n', height(resultsPosthoc));
end

fprintf('\n===== DONE =====\n');
fprintf('Output in: %s\n', outFolder);

end

function [p, test_name] = run_overall_lme_factor(T, metric, factor_name, random_factor, levels)
% Omnibus LME: metric ~ factor + (1|random). Kruskal-Wallis fallback if it does not fit.
p = NaN; test_name = '';
Tsub = T(~isnan(T.(metric)), :);
if height(Tsub) < 6; test_name = 'insufficient_n'; return; end
try
    Tlme = Tsub;
    Tlme.(factor_name)   = categorical(cellstr(Tlme.(factor_name)), levels);
    Tlme.(random_factor) = categorical(cellstr(Tlme.(random_factor)));
    lme = fitlme(Tlme, sprintf('%s ~ %s + (1|%s)', metric, factor_name, random_factor));
    anv = anova(lme);
    idx = find(strcmp(anv.Term, factor_name), 1);
    if isempty(idx); idx = height(anv); end
    p = anv.pValue(idx);
    test_name = 'LME (ParamSet effect)';
    return;
catch ME
    warning('LME (ParamSet) failed: %s. Fallback Kruskal-Wallis.', ME.message);
end
try
    p = kruskalwallis(Tsub.(metric), Tsub.(factor_name), 'off');
    test_name = 'Kruskal-Wallis (fallback)';
catch
    p = NaN; test_name = 'FAILED';
end
end

function [p, test_name, d_eff] = run_pairwise_test(xA, xB, paired)
xA = xA(~isnan(xA)); xB = xB(~isnan(xB));
p = NaN; test_name = 'insufficient_n'; d_eff = NaN;
if numel(xA) < 3 || numel(xB) < 3; return; end

nA_ok = is_normal(xA);
nB_ok = is_normal(xB);

if nA_ok && nB_ok
    try; equal_var = vartest2(xA,xB) > 0.05; catch; equal_var = false; end
    if paired
        [~,p] = ttest(xA,xB);  test_name = 'paired t-test';
    elseif equal_var
        [~,p] = ttest2(xA,xB); test_name = 't-test (equal var)';
    else
        [~,p] = ttest2(xA,xB,'Vartype','unequal'); test_name = 'Welch t-test';
    end
else
    if paired
        try;  p = signrank(xA,xB); test_name = 'Wilcoxon signed-rank';
        catch; p = ranksum(xA,xB); test_name = 'Mann-Whitney (fallback)'; end
    else
        p = ranksum(xA,xB); test_name = 'Mann-Whitney U';
    end
end

sA = std(xA,'omitnan'); sB = std(xB,'omitnan');
nA = numel(xA);         nB = numel(xB);
sp = sqrt(((nA-1)*sA^2+(nB-1)*sB^2) / max(nA+nB-2,1));
if sp > 0; d_eff = (mean(xA,'omitnan')-mean(xB,'omitnan')) / sp; end
end

function tf = is_normal(x)
x = x(~isnan(x));
if numel(x) < 4; tf = false; return; end
try;   [~,p] = swtest(x, 0.05);
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
p_adj_sorted = nan(m,1);
p_adj_sorted(si) = p_adj_v;
p_adj = nan(n,1);
p_adj(valid) = p_adj_sorted;
end

function s = pValueToStar(p)
if     p < 0.001; s = '***';
elseif p < 0.01;  s = '**';
elseif p < 0.05;  s = '*';
else;             s = 'ns';
end
end

function drawStyledBox(ax, xPos, vals, clr, bw)
if numel(vals) < 2, return; end
q1   = prctile(vals,25); q2 = prctile(vals,50); q3 = prctile(vals,75);
iqrV = q3 - q1;
wLow  = max(min(vals(vals >= q1-1.5*iqrV)), q1-1.5*iqrV);
wHigh = min(max(vals(vals <= q3+1.5*iqrV)), q3+1.5*iqrV);
outliers = vals(vals < wLow | vals > wHigh);
patch(ax, [xPos-bw/2,xPos+bw/2,xPos+bw/2,xPos-bw/2],[q1,q1,q3,q3], clr, ...
    'FaceAlpha',0.45,'EdgeColor',clr,'LineWidth',1.2);
line(ax, [xPos-bw/2,xPos+bw/2],[q2,q2], 'Color',clr*0.60,'LineWidth',2.0);
line(ax, [xPos,xPos],[wLow,q1],  'Color',clr,'LineStyle','--','LineWidth',1.0);
line(ax, [xPos,xPos],[q3,wHigh], 'Color',clr,'LineStyle','--','LineWidth',1.0);
capW = bw*0.35;
line(ax, [xPos-capW,xPos+capW],[wLow,wLow],   'Color',clr,'LineWidth',1.0);
line(ax, [xPos-capW,xPos+capW],[wHigh,wHigh], 'Color',clr,'LineWidth',1.0);
if ~isempty(outliers)
    plot(ax, repmat(xPos,numel(outliers),1), outliers, '+', ...
        'Color',clr,'MarkerSize',6,'LineWidth',1.0);
end
end

function clean = removeOutliersIQR(data, k)
data = data(isfinite(data));
if numel(data) < 4; clean = data; return; end
q1 = prctile(data,25); q3 = prctile(data,75); iqr_val = q3-q1;
clean = data(data >= q1-k*iqr_val & data <= q3+k*iqr_val);
end
