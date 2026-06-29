function NBrate_boxplot_V1vsV2_2fig()
% NBrate_boxplot_V1vsV2_2fig  NB-rate and NB-count boxplots comparing Tampere V1 vs V2.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902, NB_Tampere_DD1C_V2,
%   NB_Tampere_DD5A_V2, NB_Tampere_TUBA_V2, NB_Tampere_C10902_V2
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Tampere_DD1C_V2','NB_Tampere_DD5A_V2','NB_Tampere_TUBA_V2','NB_Tampere_C10902_V2'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% BOXPLOT: Tampere V1 vs V2 — 2 separate figures
% Figure 1 — NB rate (NBs/min)   layout 1×4 (one subplot per DIV)
% Figure 2 — Total NB count      layout 1×4 (one subplot per DIV)
% Per ogni cell line: 2 box affiancati
%   Blue   = NB_Tampere_(cellline)    [V1]
%   Orange = NB_Tampere_(cellline)_V2 [V2]
% Data point = 1 well (all MEA plates pooled).
% Inactive wells (no NB detected) → excluded + debug message.
% REQUIRED IN WORKSPACE:
%   NB_Tampere_DD1C,    NB_Tampere_DD5A,    NB_Tampere_TUBA,    NB_Tampere_C10902
%   NB_Tampere_DD1C_V2, NB_Tampere_DD5A_V2, NB_Tampere_TUBA_V2, NB_Tampere_C10902_V2

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
recordingTime_s = 600;
recordingTime_m = recordingTime_s / 60;   % 10 min → NBs/min
outlierK        = 1.5;
saveFigures     = true;
saveFolder      = fullfile(outRoot, 'all_outputs', 'NBR_V1vsV2');

if saveFigures && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

colorV1 = [0.80 0.20 0.15];  % V1 (rosso)
colorV2 = [0.00 0.65 0.25];  % V2 (verde)

%% ===== CELL LINE DEFINITIONS =====
allNames      = {'DD1C','DD5A','TUBA','C10902'};
CellLineNames = {};
for i = 1:numel(allNames)
    n  = allNames{i};
    v1 = sprintf('NB_Tampere_%s',    n);
    v2 = sprintf('NB_Tampere_%s_V2', n);
    if evalin('base',sprintf('exist(''%s'',''var'')',v1)) && ...
       evalin('base',sprintf('exist(''%s'',''var'')',v2))
        CellLineNames{end+1} = n; %#ok
    else
        if ~evalin('base',sprintf('exist(''%s'',''var'')',v1)), fprintf('[SKIP] %s not found.\n',v1); end
        if ~evalin('base',sprintf('exist(''%s'',''var'')',v2)), fprintf('[SKIP] %s not found.\n',v2); end
    end
end
if isempty(CellLineNames), error('No valid cell line with both V1 and V2.'); end

nCL  = numel(CellLineNames);
nDIV = numel(selectedDIVs);

%% ===== DATA COLLECTION =====
DATA_V1_rate  = cell(nCL,nDIV);
DATA_V2_rate  = cell(nCL,nDIV);
DATA_V1_count = cell(nCL,nDIV);
DATA_V2_count = cell(nCL,nDIV);

for c = 1:nCL
    n    = CellLineNames{c};
    NB_1 = evalin('base', sprintf('NB_Tampere_%s',    n));
    NB_2 = evalin('base', sprintf('NB_Tampere_%s_V2', n));

    for d = 1:nDIV
        DIV = selectedDIVs{d};

        cnt1 = collectNBcount(NB_1, DIV, n, 'V1', DIV);
        cnt2 = collectNBcount(NB_2, DIV, n, 'V2', DIV);

        if outlierK > 0
            cnt1 = removeOutliersIQR(cnt1, outlierK);
            cnt2 = removeOutliersIQR(cnt2, outlierK);
        end

        DATA_V1_count{c,d} = cnt1;
        DATA_V2_count{c,d} = cnt2;
        DATA_V1_rate{c,d}  = cnt1 / recordingTime_m;   % NBs/min
        DATA_V2_rate{c,d}  = cnt2 / recordingTime_m;   % NBs/min

        fprintf('[DATA] %s | %s | V1: %d active wells | V2: %d active wells\n', ...
            n, DIV, numel(cnt1), numel(cnt2));
    end
end

%% ===== FIGURE 1 — NB rate (NBs/min) =====
makeFigure(DATA_V1_rate, DATA_V2_rate, CellLineNames, selectedDIVs, ...
    colorV1, colorV2, outlierK, ...
    sprintf('NB rate (NBs/min)  [rec = %d s]', recordingTime_s), ...
    sprintf('NB Rate — Tampere V1 vs V2  (IQR k=%.1f)', outlierK), ...
    '1 data point = 1 well  |  num_NB / 10 min', ...
    'NBrate_boxplot_V1vsV2_rate', saveFigures, saveFolder);

% %% ===== FIGURE 2 — Total NB count =====
% makeFigure(DATA_V1_count, DATA_V2_count, CellLineNames, selectedDIVs, ...
%     colorV1, colorV2, outlierK, ...
%     'Total NB count', ...
%     sprintf('Total NB Count — Tampere V1 vs V2  (IQR k=%.1f)', outlierK), ...
%     '1 data point = 1 well', ...
%     'NBcount_boxplot_V1vsV2', saveFigures, saveFolder);

% LOCAL FUNCTIONS

end

function vals = collectNBcount(NB_Cell, DIV, clName, verName, divName)
% Returns NB count for each ACTIVE well.
% Wells with no NB detected → [INACTIVE] debug + excluded.
    vals = [];
    if isempty(NB_Cell) || ~isstruct(NB_Cell), return; end
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
                    clName, verName, divName, meaPlate, wellName);
                continue
            end
            vals(end+1) = nNB; %#ok
        end
    end
end

function makeFigure(DATA_1, DATA_2, CellLineNames, selectedDIVs, ...
        colorV1, colorV2, outlierK, yLabelStr, titleStr, subtitleStr, ...
        figNameBase, saveFig, saveFolder)

    nCL  = numel(CellLineNames);
    nDIV = numel(selectedDIVs);
    bw = 0.28; offV1 = -0.22; offV2 = 0.22;

    fig = figure('Color','w','Position',[60 60 350*nDIV+100 480]);
    tl  = tiledlayout(1, nDIV, 'TileSpacing','compact','Padding','compact');
%     title(tl, sprintf('%s\n%s', titleStr, subtitleStr), 'FontSize',11,'FontWeight','bold');

    for d = 1:nDIV
        ax = nexttile;
        hold(ax,'on');
        yAll = [];
        for c = 1:nCL
            d1 = DATA_1{c,d}; d2 = DATA_2{c,d};
            yAll = [yAll; d1(:); d2(:)]; %#ok
            if ~isempty(d1), drawBox(ax, c+offV1, d1, colorV1, bw); end
            if ~isempty(d2), drawBox(ax, c+offV2, d2, colorV2, bw); end
        end
        xticks(ax, 1:nCL); xticklabels(ax, CellLineNames);
        xtickangle(ax, 25); ax.XAxis.FontSize = 10;
        xlim(ax, [0.5, nCL+0.5]);
        if ~isempty(yAll) && any(isfinite(yAll))
            ylim(ax, [0, max(yAll(isfinite(yAll)))*1.15]);
        end
        ylabel(ax, yLabelStr, 'FontSize',15);
%         title(ax, selectedDIVs{d}, 'FontSize',12,'FontWeight','bold');
        grid(ax,'on'); box(ax,'on');
        if d == 1
            hV1 = patch(ax,NaN,NaN,colorV1,'EdgeColor',colorV1*0.55,'LineWidth',1.2);
            hV2 = patch(ax,NaN,NaN,colorV2,'EdgeColor',colorV2*0.55,'LineWidth',1.2);
            legend(ax,[hV1 hV2],{'V1','V2'},'Location','northeast','FontSize',10);
        end
    end

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
