function venn_fragments_fixed()
% venn_fragments_fixed  Venn diagrams of fragmented NBs (Tampere vs Doorn overlap).
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Doorn_DD1C,
%   NB_Tampere_DD5A, NB_Tampere_TUBA, NB_Tampere_C10902, NB_Doorn_DD5A,
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
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Doorn_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902','NB_Doorn_DD5A','NB_Doorn_TUBA','NB_Doorn_C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% VENN DIAGRAM — Fragment overlap: Tampere vs Doorn
% For each cell line × DIV:
%   - Overlap is computed WELL-BY-WELL (matching within the same well)
%   - Results are SUMMED across all wells before plotting
%   - One Venn diagram per cell line × DIV
% This avoids spurious cross-well matches (fragments from different wells
% at the same absolute time would never have been matched in the old
% pooled version but were at risk of it).
% REQUIRED IN WORKSPACE:
%   NB_Tampere_DD1C, NB_Doorn_DD1C  (same schema for DD5A, TUBA, C10902)
% OUTPUT FIELDS used:
%   NB.fragmentPeakTimes  (cell array, one cell per NB in the well)

%% ===== USER SETTINGS =====
selectedDIVs    = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
toleranceWindow = 0.05;   % seconds: two peaks are "same fragment" if |t1-t2| <= this
saveFigures     = true;
saveFolder      = fullfile(outRoot, 'all_outputs', 'venn');

if saveFigures && ~exist(saveFolder,'dir'), mkdir(saveFolder); end

allNames = {'DD1C','DD5A','TUBA','C10902'};

%% ===== MAIN LOOP =====
for c = 1:numel(allNames)
    n    = allNames{c};
    varT = sprintf('NB_Tampere_%s', n);
    varD = sprintf('NB_Doorn_%s',   n);

    hasT = evalin('base', sprintf('exist(''%s'',''var'')', varT));
    hasD = evalin('base', sprintf('exist(''%s'',''var'')', varD));
    if ~hasT || ~hasD
        fprintf('Skipping %s: missing %s or %s\n', n, varT, varD);
        continue
    end

    NB_T = evalin('base', varT);
    NB_D = evalin('base', varD);

    for d = 1:numel(selectedDIVs)
        DIV = selectedDIVs{d};

        % ---- well-by-well overlap, then sum ----
        [nOnly_T, nOnly_D, nShared, nWells] = ...
            computeOverlapWellByWell(NB_T, NB_D, DIV, toleranceWindow);

        if nWells == 0
            fprintf('No fragment data: %s %s\n', n, DIV);
            continue
        end

        nTotal_T = nOnly_T + nShared;
        nTotal_D = nOnly_D + nShared;
        tot      = nOnly_T + nOnly_D + nShared;
        pctS     = 0;
        if tot > 0, pctS = 100 * nShared / tot; end

        fprintf('%s | %s | wells=%d | Tampere=%d Doorn=%d Shared=%d (%.1f%%)\n', ...
            n, DIV, nWells, nTotal_T, nTotal_D, nShared, pctS);

        % ---- draw Venn ----
        fig = figure('Color','w','Position',[100 100 700 520]);
        ax  = axes(fig,'Position',[0 0 1 1]); %#ok
        axis(ax,'off'); hold(ax,'on'); axis(ax,'equal');
        xlim(ax,[-3.5 3.5]); ylim(ax,[-2.5 2.5]);

        colT = [0.80 0.20 0.15];
        colD = [0.15 0.40 0.75];   % red   — Doorn
        colO = [0.45 0.25 0.65];   % purple — overlap

        r   = 1.5;    % circle radius
        sep = 1.0;    % half-distance between centres
        cxT = -sep;  cxD = sep;  cy = 0;

        theta = linspace(0, 2*pi, 400);
        patch(ax, cxT + r*cos(theta), cy + r*sin(theta), colT, ...
            'FaceAlpha',0.30,'EdgeColor',colT,'LineWidth',2.5);
        patch(ax, cxD + r*cos(theta), cy + r*sin(theta), colD, ...
            'FaceAlpha',0.30,'EdgeColor',colD,'LineWidth',2.5);

        % numbers
        text(ax, cxT - r*0.52, 0, sprintf('%d', nOnly_T), ...
            'FontSize',22,'FontWeight','bold','Color',colT,'HorizontalAlignment','center');
        text(ax, 0, 0, sprintf('%d', nShared), ...
            'FontSize',22,'FontWeight','bold','Color',colO,'HorizontalAlignment','center');
        text(ax, cxD + r*0.52, 0, sprintf('%d', nOnly_D), ...
            'FontSize',22,'FontWeight','bold','Color',colD,'HorizontalAlignment','center');

        % labels
        text(ax, cxT - r*0.55, -r*0.82, sprintf('Tampere\n(total: %d)', nTotal_T), ...
            'FontSize',12,'Color',colT,'HorizontalAlignment','center','FontWeight','bold');
        text(ax, cxD + r*0.55, -r*0.82, sprintf('Doorn\n(total: %d)', nTotal_D), ...
            'FontSize',12,'Color',colD,'HorizontalAlignment','center','FontWeight','bold');
        % FIX: \x00B1 non viene interpretato da sprintf -> uso char(177) = ±
        text(ax, 0, r*0.78, sprintf('shared\n(%c%.2fs)', char(177), toleranceWindow), ...
            'FontSize',10,'Color',colO,'HorizontalAlignment','center');

        % title (includes well count)
%         title(ax, sprintf('%s | %s | Fragment overlap — %d wells (Tampere vs Doorn)', ...
%             n, DIV, nWells), ...
%             'FontSize',12,'FontWeight','bold','Units','normalized','Position',[0.5 0.97]);

        % annotation bottom
        annotation(fig,'textbox',[0.10 0.02 0.80 0.07], ...
            'String', sprintf('Well-by-well matching | %d/%d total unique fragments shared (%.1f%%) | tol = %.2f s', ...
            nShared, tot, pctS, toleranceWindow), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',10,'EdgeColor','none','Color',[0.3 0.3 0.3]);

        if saveFigures
            fname = sprintf('Venn_fragments_%s_%s.png', n, DIV);
            exportgraphics(fig, fullfile(saveFolder, fname), 'Resolution',200);
            fprintf('  Saved: %s\n', fname);
        end
    end
end

fprintf('\nDone.\n');

% LOCAL FUNCTIONS

end

function [totOnlyT, totOnlyD, totShared, nWells] = ...
        computeOverlapWellByWell(NB_T_Cell, NB_D_Cell, DIV, tol)
% Iterate over matching wells in both structs.
% For each well: extract fragment peaks, run greedy matching.
% Accumulate counts. Wells present in only one detector are included
% (the other contributes 0 peaks → all fragments count as exclusive).

totOnlyT = 0;  totOnlyD = 0;  totShared = 0;  nWells = 0;

% FIX: union su cell array Nx2 con 'rows' è instabile.
% Codifico (plate, well) come singola stringa "plate||well", faccio
% unique sul vettore di stringhe, poi spacchetto in cell array Nx2.
kT = collectWellKeys(NB_T_Cell, DIV);
kD = collectWellKeys(NB_D_Cell, DIV);

encT = {};
if ~isempty(kT), encT = strcat(kT(:,1), '||', kT(:,2)); end
encD = {};
if ~isempty(kD), encD = strcat(kD(:,1), '||', kD(:,2)); end

encAll = unique([encT; encD]);
if isempty(encAll), return; end

keys = cell(numel(encAll), 2);
for i = 1:numel(encAll)
    parts = strsplit(encAll{i}, '||');
    keys(i,:) = parts;
end

for k = 1:size(keys,1)
    mp  = keys{k,1};
    wl  = keys{k,2};

    peaksT = getWellPeaks(NB_T_Cell, mp, DIV, wl);
    peaksD = getWellPeaks(NB_D_Cell, mp, DIV, wl);

    if isempty(peaksT) && isempty(peaksD), continue; end
    nWells = nWells + 1;

    [oT, oD, sh] = greedyMatch(peaksT, peaksD, tol);
    totOnlyT  = totOnlyT  + oT;
    totOnlyD  = totOnlyD  + oD;
    totShared = totShared + sh;
end
end

function keys = collectWellKeys(NB_Cell, DIV)
% Returns Nx2 cell array of {plate, well} for all wells at this DIV
% FIX: inizializzo come cell(0,2) invece di {} per coerenza dimensionale
keys = cell(0,2);
if isempty(NB_Cell) || ~isstruct(NB_Cell), return; end
plates = fieldnames(NB_Cell);
for p = 1:numel(plates)
    mp = plates{p};
    if ~isfield(NB_Cell.(mp), DIV), continue; end
    wells = fieldnames(NB_Cell.(mp).(DIV));
    for w = 1:numel(wells)
        keys(end+1,:) = {mp, wells{w}}; %#ok
    end
end
end

function peaks = getWellPeaks(NB_Cell, mp, DIV, wl)
% Extract and concatenate all fragmentPeakTimes for one well
peaks = [];
if isempty(NB_Cell) || ~isstruct(NB_Cell), return; end
if ~isfield(NB_Cell, mp),                  return; end
if ~isfield(NB_Cell.(mp), DIV),            return; end
if ~isfield(NB_Cell.(mp).(DIV), wl),       return; end

NB = NB_Cell.(mp).(DIV).(wl);
if isempty(NB) || ~isstruct(NB),                           return; end
if ~isfield(NB,'fragmentPeakTimes') || ...
   isempty(NB.fragmentPeakTimes),                          return; end

% pre-conto per pre-allocare
nTot = 0;
for i = 1:numel(NB.fragmentPeakTimes)
    nTot = nTot + numel(NB.fragmentPeakTimes{i});
end
peaks = zeros(nTot,1);
idx = 0;
for i = 1:numel(NB.fragmentPeakTimes)
    fp = NB.fragmentPeakTimes{i};
    if ~isempty(fp)
        m = numel(fp);
        peaks(idx+1:idx+m) = fp(:);
        idx = idx + m;
    end
end
peaks = sort(peaks(1:idx));
end

function [nOnlyA, nOnlyB, nShared] = greedyMatch(peaksA, peaksB, tol)
% Greedy 1-to-1 matching within tolerance window.
% FIX: gestisco i casi vuoti fuori dal loop per chiarezza.

if isempty(peaksA) && isempty(peaksB)
    nOnlyA = 0; nOnlyB = 0; nShared = 0; return
end
if isempty(peaksB)
    nOnlyA  = numel(peaksA);
    nOnlyB  = 0;
    nShared = 0;
    return
end
if isempty(peaksA)
    nOnlyA  = 0;
    nOnlyB  = numel(peaksB);
    nShared = 0;
    return
end

matchedA = false(size(peaksA));
matchedB = false(size(peaksB));

for i = 1:numel(peaksA)
    diffs = abs(peaksB - peaksA(i));
    [minD, idx] = min(diffs);
    if minD <= tol && ~matchedB(idx)
        matchedA(i)   = true;
        matchedB(idx) = true;
    end
end

nShared = sum(matchedA);
nOnlyA  = sum(~matchedA);
nOnlyB  = sum(~matchedB);
end
