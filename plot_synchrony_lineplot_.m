function plot_synchrony_lineplot_()
% plot_synchrony_lineplot_  Line plot (95% CI) of mean active channels per NB across DIVs.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902, NB_Tampere_DD1C, NB_Tampere_DD5A,
%   NB_Tampere_TUBA, NB_Tampere_C10902
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902','NB_Tampere_DD1C','NB_Tampere_DD5A','NB_Tampere_TUBA','NB_Tampere_C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% PLOT_SYNCHRONY_LINEPLOT.m
% Line plot con 95% CI: numero medio di canali attivi per NB
% per ogni cell line attraverso i DIV selezionati.
% METRICA: per ogni well -> media canali unici in NB.CHANNELSFIN
%          -> per DIV: media ± 95% CI tra well (dopo IQR removal)
% LABEL:   valore medio arrotondato (corrisponde al punto Y)
% REQUIRED IN WORKSPACE:
% NB_Tampere_DD1C, NB_Tampere_DD5A, NB_Tampere_TUBA, NB_Tampere_C10902

%% USER SETTINGS
selectedDIVs  = {'DIV53','DIV60','DIV67','DIV83','DIV88'};
divLabels     = [53, 60, 67, 83, 90];
cellLineNames = {'DD1C','DD5A','TUBA','C10902'};

outlierK = 1.5;

saveFigures = true;
saveFolder  = fullfile(outRoot, 'all_outputs', 'synchrony');
if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

CL_colors = [0.00 0.45 0.74;  % DD1C - blu
             0.85 0.33 0.10;  % DD5A - arancio
             0.47 0.67 0.19;  % TUBA - verde
             0.49 0.18 0.56]; % C10902 - viola

%% LOAD STRUCTS
NBstructs = {NB_Tampere_DD1C, NB_Tampere_DD5A, NB_Tampere_TUBA, NB_Tampere_C10902};

nCL  = numel(cellLineNames);
nDIV = numel(selectedDIVs);

%% ESTRAI DATI
allChanPerNB = cell(nCL, nDIV);

for c = 1:nCL
    NBstruct   = NBstructs{c};
    plateNames = fieldnames(NBstruct);

    for d = 1:nDIV
        DIV  = selectedDIVs{d};
        vals = [];

        for p = 1:numel(plateNames)
            mea = plateNames{p};
            if ~isfield(NBstruct.(mea), DIV), continue; end

            wellNames = fieldnames(NBstruct.(mea).(DIV));
            for w = 1:numel(wellNames)
                NB = NBstruct.(mea).(DIV).(wellNames{w});

                if ~isfield(NB,'starts') || isempty(NB.starts) || all(isnan(NB.starts))
                    continue
                end
                if ~isfield(NB,'CHANNELSFIN') || isempty(NB.CHANNELSFIN)
                    continue
                end

                nNB = numel(NB.CHANNELSFIN);
                ch  = zeros(1, nNB);
                for k = 1:nNB
                    if ~isempty(NB.CHANNELSFIN{k})
                        ch(k) = height(NB.CHANNELSFIN{k});
                    end
                end
                ch = ch(ch > 0);
                if isempty(ch), continue; end

                vals(end+1) = mean(ch); %#ok
            end
        end

        allChanPerNB{c,d} = removeOutliersIQR(vals(:), outlierK);
    end
end

%% CALCOLA MEDIA E 95% CI per ogni cell line x DIV
meanVals = NaN(nCL, nDIV);
ciVals   = NaN(nCL, nDIV);
nVals    = zeros(nCL, nDIV);

for c = 1:nCL
    for d = 1:nDIV
        v = allChanPerNB{c,d};
        n = numel(v);
        nVals(c,d) = n;
        if n >= 1
            meanVals(c,d) = mean(v);
        end
        if n >= 2
            t_crit      = tinv(0.975, n-1);
            ciVals(c,d) = t_crit * std(v) / sqrt(n);
        end
    end
end

%% FIGURA
fig = figure('Color','w','Position',[100 100 680 520], ...
    'Name','Network Synchrony - Line Plot');

ax = axes(fig);
hold(ax,'on');

hLines = gobjects(nCL,1);

for c = 1:nCL
    clr = CL_colors(c,:);
    

    valid = ~isnan(meanVals(c,:));
    xv = divLabels(valid);
    yv = meanVals(c, valid);
    sv = ciVals(c, valid);

    if isempty(xv), continue; end

    % Banda 95% CI shaded
    x_patch = [xv, fliplr(xv)];
    y_patch = [yv + sv, fliplr(yv - sv)];
    y_patch(isnan(y_patch)) = [yv(isnan(sv)), fliplr(yv(isnan(sv)))];

    fill(ax, x_patch, y_patch, clr, ...
        'FaceAlpha', 0.18, 'EdgeColor','none');

    % Linea principale
    hLines(c) = plot(ax, xv, yv, '-o', ...
        'Color',           clr, ...
        'LineWidth',       1.8, ...
         ...
        'MarkerSize',      6, ...
        'MarkerFaceColor', clr )
        

end

% Label: valore medio (corrisponde al punto Y), ordine DD1C>DD5A>TUBA>C10902
allVals_flat = meanVals(:) + ciVals(:);
y_top  = max(allVals_flat(~isnan(allVals_flat)));
y_top  = y_top + range_y_offset(meanVals) * 0.10;
lineH  = range_y_offset(meanVals) * 0.09;
nOrder = [1, 2, 3, 4];

for d = 1:nDIV
    xpos = divLabels(d);
    for row = 1:nCL
        c    = nOrder(row);
        ypos = y_top + (nCL - row) * lineH;

        if ~isnan(meanVals(c,d))
            text(ax, xpos, ypos, ...
                sprintf('m=%.1f', meanVals(c,d)), ...   % valore medio = punto Y
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   'middle', ...
                'FontSize',            7, ...
                'Color',               CL_colors(c,:) * 0.75);
        end
    end
end

% Espandi asse Y per fare spazio alle label
y_min = min(meanVals(:) - ciVals(:), [], 'all', 'omitnan');
y_max = y_top + nCL * lineH + range_y_offset(meanVals) * 0.05;
ylim(ax, [y_min * 0.92, y_max]);

% Assi e stile
set(ax, 'XTick', divLabels, ...
    'XTickLabel', arrayfun(@(x) sprintf('DIV%d',x), divLabels, 'UniformOutput',false), ...
    'FontSize', 10, 'Box','on', 'GridAlpha',0.25);
grid(ax,'on');
xlabel(ax, 'DIV',                                     'FontSize', 15, 'FontWeight','bold');
ylabel(ax, 'Active channels per NB (mean ± 95% CI)', 'FontSize', 15, 'FontWeight','bold');
% title(ax,  'Network Synchrony — Active Channels per NB', ...
%     'FontSize', 12, 'FontWeight','bold');

validCL = find(any(~isnan(meanVals), 2));
legend(ax, hLines(validCL), cellLineNames(validCL), ...
    'Location','best', 'FontSize', 9, 'Box','on');

hold(ax,'off');

%% SALVATAGGIO
if saveFigures
    fname = 'fig_synchrony_channels_per_NB_lineplot_D88.png';
    exportgraphics(fig, fullfile(saveFolder, fname), 'Resolution', 300);
    fprintf('Figura salvata: %s\n', fname);
end

fprintf('\n--- Riepilogo (mean ± 95%% CI, IQR k=%.1f) ---\n', outlierK);
for c = 1:nCL
    fprintf('\n%s:\n', cellLineNames{c});
    for d = 1:nDIV
        fprintf('  %s: mean=%.2f | 95%%CI=%.2f | n_wells=%d\n', ...
            selectedDIVs{d}, meanVals(c,d), ciVals(c,d), nVals(c,d));
    end
end

%% LOCAL FUNCTIONS
end

function clean = removeOutliersIQR(data, k)
    data = data(isfinite(data));
    if numel(data) < 4, clean = data; return; end
    q1      = prctile(data, 25);
    q3      = prctile(data, 75);
    iqr_val = q3 - q1;
    clean   = data(data >= q1 - k*iqr_val & data <= q3 + k*iqr_val);
end

function r = range_y_offset(M)
    v = M(~isnan(M));
    if isempty(v) || range(v)==0
        r = 1;
    else
        r = range(v);
    end
end
