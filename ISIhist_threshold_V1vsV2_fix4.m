function ISIhist_threshold_V1vsV2_fix4()
% ISIhist_threshold_V1vsV2_fix4  Log ISI-N histograms per well with Tampere V1 vs V2 thresholds.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, Params
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','Params'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% ISI-N LOG HISTOGRAM per well — Threshold V1 (red) vs V2 (green)
% Per ogni cell line × DIV selezionato × well disponibile:
%   - costruisce il log ISI-N histogram con construct_ISIhist
%   - calcola la soglia con logica V1 (rosso): primo local min dopo il max
%   - calcola la soglia con logica V2 (verde): min tra i primi 2 peak rilevanti
%   - plotta l'istogramma + le 2 linee verticali
% Layout: una figura per cell line × DIV, subplot 4×4 (one per well)
% REQUIRED IN WORKSPACE:
%   DD1C.meaPlate.DIV.well   (table con colonne Channel, Time)
%   DD5A, TUBA, C10902 — stesso schema
% REQUIRED ON PATH:
%   construct_ISIhist.m

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};

% V2 peak filter: picchi < minPeakFrac * max(HH) vengono ignorati
minPeakFrac = 0.10;

nPerFig     = 10;    % max subplot per figura
saveFigures = true;
saveFolder  = fullfile(outRoot, 'all_outputs', 'ISIhist_V1vsV2');

if saveFigures && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

%% ===== CELL LINE DEFINITIONS =====
allNames = {'DD1C','DD5A','TUBA','C10902'};
CellLineNames = {};
for i = 1:numel(allNames)
    if evalin('base', sprintf('exist(''%s'',''var'')', allNames{i}))
        CellLineNames{end+1} = allNames{i}; %#ok<AGROW>
    end
end
if isempty(CellLineNames), error('Nessuna struttura dati trovata nel workspace.'); end

%% ===== LOOP CELL LINE × DIV =====
for c = 1:numel(CellLineNames)
    clName = CellLineNames{c};
    DATA   = evalin('base', clName);

    for d = 1:numel(selectedDIVs)
        DIV = selectedDIVs{d};

        % raccogli tutti i well disponibili per questo DIV
        WellList = {};  % {meaPlate, wellName, spikeTable}
        plates = fieldnames(DATA);
        for p = 1:numel(plates)
            mp = plates{p};
            if ~isfield(DATA.(mp), DIV), continue; end
            wells = fieldnames(DATA.(mp).(DIV));
            for w = 1:numel(wells)
                wn  = wells{w};
                tbl = DATA.(mp).(DIV).(wn);
                if isempty(tbl) || height(tbl) < Params.LogISINHist.N + 1, continue; end
                if ~ismember('Time',tbl.Properties.VariableNames), continue; end
                WellList{end+1,1} = mp;  %#ok<AGROW>
                WellList{end,2}   = wn;
                WellList{end,3}   = tbl;
            end
        end

        if isempty(WellList)
            fprintf('No data: %s %s\n', clName, DIV);
            continue
        end

        nWells   = size(WellList,1);
        nFigures = ceil(nWells / nPerFig);

        for f = 1:nFigures
            idxStart = (f-1)*nPerFig + 1;
            idxEnd   = min(f*nPerFig, nWells);

            fig = figure('Color','w','Position',[50 50 1700 1100]);
            tl  = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
%             title(tl, sprintf('%s  |  %s  |  page %d/%d', clName, DIV, f, nFigures), ...
%                 'FontSize',12,'FontWeight','bold');

            for k = idxStart:idxEnd
                mp  = WellList{k,1};
                wn  = WellList{k,2};
                tbl = WellList{k,3};

                ax = nexttile;
                hold(ax,'on');

                SpikeTimes = tbl.Time;

                % --- costruisci ISI histogram ---
                [Steps, ~, HH] = construct_ISIhist(SpikeTimes, ...
                    Params.LogISINHist.N, Params.LogISINHist.Span, ...
                    Params.LogISINHist.BinsPerDec, Params.LogISINHist.Smoother);

                if isempty(Steps) || isempty(HH)
%                     title(ax, sprintf('%s\\%s — empty hist', mp, wn), ...
%                         'Interpreter','none','FontSize',8);
                    axis(ax,'off'); continue
                end

                HHrnd = round(HH, 9);
                xIdx  = 1:numel(HHrnd);

                % plot histogram — patch stairs-style (compatibile con log scale)
                % costruisce i vertici del profilo a gradini
                nBins  = numel(Steps);
                % edges sinistro/destro di ogni bin in log space
                if nBins > 1
                    logS   = log10(Steps);
                    dlog   = diff(logS);
                    edges  = [logS(1)-dlog(1)/2, logS(1:end-1)+dlog/2, logS(end)+dlog(end)/2];
                    edges  = 10.^edges;
                else
                    edges = [Steps*0.5, Steps*2];
                end
                % costruisce x,y per patch a gradini
                xPatch = zeros(1, 4*nBins);
                yPatch = zeros(1, 4*nBins);
                for ib = 1:nBins
                    base = (ib-1)*4;
                    xPatch(base+1) = edges(ib);
                    xPatch(base+2) = edges(ib);
                    xPatch(base+3) = edges(ib+1);
                    xPatch(base+4) = edges(ib+1);
                    yPatch(base+1) = 0;
                    yPatch(base+2) = HHrnd(ib);
                    yPatch(base+3) = HHrnd(ib);
                    yPatch(base+4) = 0;
                end
                patch(ax, xPatch, yPatch, [0.75 0.85 0.95], ...
                    'EdgeColor',[0.45 0.60 0.80],'LineWidth',0.5,'FaceAlpha',0.9);
                set(ax,'XScale','log');

                % ---- THRESHOLD V1: primo local min dopo il max ----
                thrV1_ms = NaN;
                [~, MaxInd] = max(HHrnd);
                Localmins   = find(islocalmin(HHrnd));
                Localmins   = Localmins(Localmins > MaxInd);
                if ~isempty(Localmins)
                    ThreshIdx_V1 = Localmins(1);
                    thrV1_ms     = Steps(ThreshIdx_V1);
                end

                % ---- THRESHOLD V2: min assoluto nella valley tra i 2 peak principali ----
                % usa findpeaks con filtro su altezza E prominenza per evitare
                % picchi spuri da smoothing noise
                thrV2_ms     = NaN;
                minPeakH     = minPeakFrac * max(HHrnd);
                minProminence = 0.15 * max(HHrnd);
                [~, peakIdx] = findpeaks(double(HHrnd), ...
                    'MinPeakHeight',    minPeakH, ...
                    'MinPeakProminence',minProminence);
                if numel(peakIdx) >= 2
                    p1 = peakIdx(1);  p2 = peakIdx(2);
                    if p2 > p1 + 1
                        % minimo GLOBALE ASSOLUTO nella valley (esclusi i peak)
                        segment      = HHrnd(p1+1 : p2-1);
                        [~, relIdx]  = min(segment);
                        ThreshIdx_V2 = p1 + relIdx;
                        thrV2_ms     = Steps(ThreshIdx_V2);
                    end
                end

                % ---- plot threshold lines ----
                ylims = ylim(ax);
                if ~isnan(thrV1_ms)
                    xline(ax, thrV1_ms, 'r-', 'LineWidth', 1.8, ...
                        'Label', sprintf('V1: %.1fms', thrV1_ms*1000), ...
                        'LabelVerticalAlignment','bottom', ...
                        'FontSize',7,'HandleVisibility','off');
                end
                if ~isnan(thrV2_ms)
                    xline(ax, thrV2_ms, 'g-', 'LineWidth', 1.8, ...
                        'Label', sprintf('V2: %.1fms', thrV2_ms*1000), ...
                        'LabelVerticalAlignment','top', ...
                        'FontSize',7,'HandleVisibility','off');
                end

                % ---- legenda solo primo subplot ----
                if k == idxStart
                    hV1 = xline(ax, NaN, 'r-', 'LineWidth',1.8);
                    hV2 = xline(ax, NaN, 'g-', 'LineWidth',1.8);
                    legend(ax,[hV1 hV2],{'Thr V1','Thr V2'}, ...
                        'Location','northeast','FontSize',8);
                end

                xlabel(ax,'ISI-N (s)','FontSize',15);
                ylabel(ax,'Norm. count','FontSize',15);
%                 title(ax, sprintf('%s \\ %s', mp, wn), ...
%                     'Interpreter','none','FontSize',8,'FontWeight','bold');
                grid(ax,'on'); box(ax,'on');
            end

            % spegni tile vuoti
            nUsed = idxEnd - idxStart + 1;
            for z = nUsed+1:nPerFig
                ax = nexttile; axis(ax,'off');
            end

            if saveFigures
                fname = sprintf('ISIhist_%s_%s_page%02d.png', clName, DIV, f);
                exportgraphics(fig, fullfile(saveFolder,fname),'Resolution',200);
                fprintf('Saved: %s\n', fname);
            end
        end % figure pages
    end % DIV
end % cell line

fprintf('Done.\n');

end
