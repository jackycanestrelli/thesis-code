function radar_S_inner_pentagon()
% radar_S_inner_pentagon  Radar / spider chart of per-cell-line activity metrics (from metrics CSV).
%
% INPUT : reads metric CSV/MAT files; output paths configurable via RESULTS_ROOT
%         in the base workspace (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% RADAR SPIDER CHART
%  USER SETTINGS
%selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
% Root holding the per-parameter-set metric CSVs (paramsX/metrics_per_well.csv).
% Override by setting METRICS_ROOT in the base workspace if the drive changes.
if evalin('base', "exist('METRICS_ROOT','var')")
    metricsRoot = evalin('base','METRICS_ROOT');
else
    metricsRoot = 'S:\met_narkilahti_neuro_sto-3700\MEA_data_internship\Jacopo\TEST SCRIPT APP';
end
csvPath = fullfile(metricsRoot, 'params2', 'metrics_per_well.csv');

% --- Cell lines -----------------------------------------------------------
% celllines : exact strings as they appear in the CSV  CellLine column
% cllabels  : display names shown in the legend (free text, one per cell line)
celllines = {'DD1C',  'DD5A',  'TUBA',    'C10902'};
cllabels  = {'DD1C (severe)', 'DD5A (mild)', 'TUBA (ctrl)', 'C10902 (ctrl)'};

% --- DIV to plot (ONE value) ----------------------------------------------
plotDIVs = [53, 60, 67, 83, 90];   % one radar per DIV

% --- Detector for NB-derived metrics (NBR, FragPrev, NBD, CVIBI) ----------
detectorName = 'TampereV1';

% --- Metrics --------------------------------------------------------------
% metricNames  : must match CSV column names exactly
% metricLabels : display names on the radar axes
metricNames  = {'MFR',  'SBR',     'NBR',     'FragPrev',    'NBD',   'CV_IBI'};
metricLabels = {'MFR',  'SEB rate', 'NB Rate', 'Frag. Prev.', 'NB Dur','CV IBI'};

% --- Per-metric scale (raw units) -----------------------------------------
% scaleMin : value mapped to the CENTRE of the radar (innermost ring)
% scaleMax : value mapped to the OUTER ring
%            MFR   SBR   NBR   FragPrev  NBD   CVIBI
scaleMin = [  0,    0,    10,    0,        0,    0  ];
scaleMax = [  20,   105,    0,    1,        30,   2  ];

% --- Number of grid rings (pentagon layers) --------------------------------
nRings = 5;   % includes the inner pentagon (scaleMin) and outer ring (scaleMax)

% --- Inner pentagon radius -------------------------------------------------
% Fixed visual radius (0..1) where scaleMin is placed for every metric.
% A cell line whose value equals scaleMin sits on this inner pentagon,
% NOT at the centre. Recommended: 0.10 .. 0.20
innerPentagonRadius = 0.15;

% --- IQR outlier removal --------------------------------------------------
enableOutlierRemoval = false;
outlierK             = 1.5;

% --- Visual ---------------------------------------------------------------
palette = [
    0.70  0.10  0.10;   % DD1C  - red
    0.95  0.55  0.10;   % DD5A  - orange
    0.10  0.30  0.70;   % TUBA  - blue
    0.30  0.60  0.85;   % C10902- light blue
];
alphaFill  = 0.15;

% --- Save -----------------------------------------------------------------
saveFigure = true;
saveFolder = fullfile(outRoot, 'all_outputs', 'radar');   % figures saved locally under RESULTS_ROOT
if saveFigure && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

%  LOAD CSV

fprintf('Loading %s ...\n', csvPath);
Tall = readtable(csvPath, 'TextType','string');
Tall.CellLine = string(Tall.CellLine);
Tall.Detector = string(Tall.Detector);
Tall.Well     = string(Tall.Well);
Tall.Plate    = string(Tall.Plate);

% One radar per DIV; the cell-line legend is shown only on the central radar.
nDIVrad = numel(plotDIVs);
for di = 1:nDIVrad
    plotDIV    = plotDIVs(di);
    showLegend = (di == ceil(nDIVrad/2));

    T = Tall(Tall.DIV == plotDIV, :);
    if isempty(T)
        warning('No data found for DIV %d. Skipping.', plotDIV);
        continue
    end
    fprintf('DIV %d selected: %d rows\n', plotDIV, height(T));

    % Separate detector-specific sub-table for NB metrics
    Tnb  = T(T.Detector == string(detectorName), :);
    Tmfr = T;   % MFR / SBR / SBD use all rows regardless of detector

nKept = height(Tnb);
nAll  = height(T);
fprintf('Detector %s: %d / %d rows kept for NB metrics (%.1f%%).\n', ...
        detectorName, nKept, nAll, 100*nKept/nAll);

% Check columns
missingCols = setdiff(metricNames, T.Properties.VariableNames);
if ~isempty(missingCols)
    error('Missing columns in CSV: %s.\nRegenerate the CSV.', strjoin(missingCols,', '));
end

%  AGGREGATE: mean across wells per cell line (for the selected DIV)
nCL = numel(celllines);
nM  = numel(metricNames);

mfrMetrics = {'MFR','SBR','SBD'};
rawAvg = nan(nCL, nM);   % [cellLine x metric]

for ic = 1:nCL
    cl = celllines{ic};
    for im = 1:nM
        mname = metricNames{im};
        if ismember(mname, mfrMetrics)
            sub = Tmfr(Tmfr.CellLine == string(cl), :);
            % deduplicate wells (take first occurrence per well)
            [~, ia] = unique(sub.Well);
            sub = sub(ia, :);
        else
            sub = Tnb(Tnb.CellLine == string(cl), :);
        end
        if isempty(sub) || ~ismember(mname, sub.Properties.VariableNames)
            continue
        end
        vals = sub.(mname);
        vals = vals(~isnan(vals));
        if isempty(vals), continue; end
        if enableOutlierRemoval
            vals = removeOutliersIQR(vals, outlierK);
            if isempty(vals), continue; end
        end
        rawAvg(ic, im) = mean(vals, 'omitnan');
    end
end

%  PRINT SUMMARY TABLE (raw values)
fprintf('\n%-12s', 'CellLine');
for im = 1:nM, fprintf('  %-10s', metricLabels{im}); end
fprintf('\n%s\n', repmat('-', 1, 12 + nM*12));
for ic = 1:nCL
    fprintf('%-12s', cllabels{ic});
    for im = 1:nM
        if isnan(rawAvg(ic,im))
            fprintf('  %-10s', 'NaN');
        else
            fprintf('  %-10.3f', rawAvg(ic,im));
        end
    end
    fprintf('\n');
end

%  PLOT
mainTitle = sprintf('DIV %d', plotDIV);   % header at the top of each radar
if enableOutlierRemoval
    outlierStr = sprintf('IQR k=%.1f', outlierK);
else
    outlierStr = 'no outlier removal';
end
subTitle  = '';

fig = drawRawRadar(rawAvg, metricLabels, cllabels, palette, alphaFill, ...
                   mainTitle, subTitle, scaleMin, scaleMax, nRings, innerPentagonRadius, showLegend);

if saveFigure
    fname = sprintf('RadarRAW_DIV%d_%s', plotDIV, detectorName);
    exportgraphics(fig, fullfile(saveFolder, [fname '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(saveFolder, [fname '.pdf']), 'ContentType', 'vector');
    fprintf('Saved %s\n', fname);
end
end   % per-DIV radar loop
fprintf('Done.\n');

%  LOCAL FUNCTIONS

end

function fig = drawRawRadar(dataMat, metricLabels, cllabels, palette, alphaFill, ...
                             mainTitle, subTitle, scaleMin, scaleMax, nRings, rInner, showLegend)
if nargin < 12 || isempty(showLegend), showLegend = true; end
% Draw a radar chart with raw (un-normalised) values, using an INNER PENTAGON
% as the baseline (scaleMin) instead of a single centre point.
% Layout:
%   - Inner pentagon at fixed visual radius rInner  (= where scaleMin sits)
%   - Outer ring at visual radius 1.0              (= where scaleMax sits)
%   - Intermediate rings between rInner and 1.0, spaced linearly.
%   - Each metric has its own scaleMin / scaleMax (per-spoke), but they all
%     share the same visual rInner..1.0 envelope.
% Value-to-radius mapping:
%   r_disp = rInner + (1 - rInner) * (value - scaleMin) / (scaleMax - scaleMin)
%   then clipped to [rInner, 1].
% NaN values are mapped to rInner (sit on the inner pentagon).

nCL = size(dataMat, 1);
nM  = size(dataMat, 2);
if size(palette,1) < nCL, palette = lines(nCL); end

% Sanity: rInner must be in (0, 1)
if rInner <= 0 || rInner >= 1
    warning('innerPentagonRadius must be in (0,1). Using 0.15.');
    rInner = 0.15;
end

% Angles: start from top (-pi/2), go clockwise
angles = linspace(0, 2*pi, nM+1);
angles = angles(1:nM) - pi/2;

fig = figure('Color','w','Position',[100 100 920 920]);
ax  = axes('Parent',fig,'Visible','off','Position',[0.05 0.10 0.90 0.78]);
hold(ax,'on');
axis(ax,'equal');

% ----- Grid rings ---------------------------------------------------------
% nRings rings between rInner (innermost) and 1.0 (outermost), spaced linearly
ringRadii = linspace(rInner, 1.0, nRings);

for ig = 1:nRings
    r  = ringRadii(ig);
    th = [angles, angles(1)];
    [xi, yi] = pol2cart(th, r * ones(1, nM+1));
    if ig == nRings
        % outer ring (scaleMax): bold black
        plot(ax, xi, yi, 'Color',[0.10 0.10 0.10], 'LineWidth', 1.6);
    elseif ig == 1
        % inner pentagon (scaleMin): dashed dark
        plot(ax, xi, yi, '--', 'Color',[0.35 0.35 0.35], 'LineWidth', 1.2);
    else
        plot(ax, xi, yi, 'Color',[0.78 0.78 0.78], 'LineWidth', 0.8);
    end
end

% ----- Spokes (drawn from inner pentagon to outer ring only) --------------
% This emphasises that the inner pentagon is the actual baseline.
for im = 1:nM
    [xi0, yi0] = pol2cart(angles(im), rInner);
    [xi1, yi1] = pol2cart(angles(im), 1.0);
    plot(ax, [xi0 xi1], [yi0 yi1], 'Color',[0.68 0.68 0.68], 'LineWidth', 1);
end

% ----- Per-spoke tick labels (scaleMin near inner, scaleMax near outer) ---
% scaleMin label: placed just INSIDE the inner pentagon, on each spoke
% scaleMax label: placed just OUTSIDE the outer ring, on each spoke
offsetIn  = 0.08;   % how far inside the inner pentagon to place the min label
offsetOut = 0.06;   % how far outside the outer ring to place the max label

for im = 1:nM
    a = angles(im);

    % --- scaleMin label (inside the inner pentagon) ---
    rMinLabel = max(0.02, rInner - offsetIn);
    [xlMin, ylMin] = pol2cart(a, rMinLabel);
    text(ax, xlMin, ylMin, sprintf('%.3g', scaleMin(im)), ...
        'FontSize', 12, 'Color',[0.30 0.30 0.30], ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'BackgroundColor','w','Margin',1);

    % --- scaleMax label (outside the outer ring) ---
    rMaxLabel = 1.0 + offsetOut;
    [xlMax, ylMax] = pol2cart(a, rMaxLabel);
    text(ax, xlMax, ylMax, sprintf('%.3g', scaleMax(im)), ...
        'FontSize', 12, 'Color',[0.10 0.10 0.10], 'FontWeight','bold', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'BackgroundColor','w','Margin',1);
end

% ----- Intermediate ring tick (along first spoke only, for reference) -----
% Show the intermediate values along the first spoke so the user has a
% sense of the linear scale between min and max.
for ig = 2:(nRings-1)
    r    = ringRadii(ig);
    frac = (r - rInner) / (1 - rInner);    % 0..1 along min..max
    val1 = scaleMin(1) + frac * (scaleMax(1) - scaleMin(1));
    [xl, yl] = pol2cart(angles(1) - 0.18, r);
    text(ax, xl, yl, sprintf('%.3g', val1), ...
        'FontSize', 11, 'Color',[0.55 0.55 0.55], ...
        'HorizontalAlignment','center','VerticalAlignment','middle');
end

% ----- Axis (metric) labels ----------------------------------------------
labelR = 1.40;
for im = 1:nM
    [xi, yi] = pol2cart(angles(im), labelR);
    % Align the text away from the centre on each side so long labels
    % (e.g. "Network synchrony") grow outward instead of into the figure edge.
    ca = cos(angles(im));
    if     ca >  0.30, hAlign = 'left';
    elseif ca < -0.30, hAlign = 'right';
    else,              hAlign = 'center';
    end
    text(ax, xi, yi, metricLabels{im}, ...
        'HorizontalAlignment',hAlign,'VerticalAlignment','middle', ...
        'FontSize', 15, 'FontWeight','bold', 'Color',[0.15 0.15 0.15]);
end

% ----- Polygons -----------------------------------------------------------
handles = gobjects(nCL, 1);
for ic = 1:nCL
    row   = dataMat(ic, :);
    rdisp = nan(1, nM);
    for im = 1:nM
        v = row(im);
        if isnan(v)
            rdisp(im) = rInner;   % missing -> inner pentagon (baseline)
        else
            rng = scaleMax(im) - scaleMin(im);
            if rng == 0, rng = 1; end
            frac = (v - scaleMin(im)) / rng;
            frac = max(0, min(1, frac));
            rdisp(im) = rInner + (1 - rInner) * frac;
        end
    end
    th = [angles, angles(1)];
    rr = [rdisp, rdisp(1)];
    [xi, yi] = pol2cart(th, rr);
    fill(ax, xi, yi, palette(ic,:), 'FaceAlpha', alphaFill, 'EdgeColor','none');
    handles(ic) = plot(ax, xi, yi, '-o', ...
        'Color', palette(ic,:), 'LineWidth', 2.2, ...
        'MarkerSize', 6, 'MarkerFaceColor', palette(ic,:), ...
        'DisplayName', cllabels{ic});
end

% ----- Title & subtitle ---------------------------------------------------
annotation(fig,'textbox',[0 0.93 1 0.06], ...
    'String', mainTitle, ...
    'FontSize',18,'FontWeight','bold', ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FitBoxToText','off');
annotation(fig,'textbox',[0 0.88 1 0.05], ...
    'String', subTitle, ...
    'FontSize',13,'FontWeight','normal','Color',[0.40 0.40 0.40], ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'EdgeColor','none','FitBoxToText','off');

% ----- Legend (only on the central radar) ---------------------------------
if showLegend
    legend(handles, cllabels, 'Location','southoutside', ...
           'Orientation','horizontal','FontSize',15,'Box','off');
end

xlim(ax, [-1.85 1.85]);
ylim(ax, [-1.85 1.85]);
end

function clean = removeOutliersIQR(data, k)
data = data(~isnan(data));
if numel(data) < 4, clean = data; return; end
q1  = prctile(data, 25);
q3  = prctile(data, 75);
iqr = q3 - q1;
clean = data(data >= q1 - k*iqr & data <= q3 + k*iqr);
end
