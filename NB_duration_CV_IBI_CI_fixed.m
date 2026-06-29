function NB_duration_CV_IBI_CI_fixed()
% NB_duration_CV_IBI_CI_fixed  Mean NB duration line plot and CV/LvR inter-burst-interval boxplots.
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

% MEAN NETWORK BURST DURATION — line plot (4 CL x 4 DIV)
% CV INTER-BURST INTERVAL — boxplot (4 DIV, 4 CL per panel)
% LvR INTER-BURST INTERVAL — boxplot (4 DIV, 4 CL per panel)
% Formula LvR (Shinomoto et al.):
%   LvR = 3/(n-1) * sum_{i=1}^{n-1} (1 - 4*Ii*Ii+1/(Ii+Ii+1)^2) * (1 + 4R/(Ii+Ii+1))
% Usa strutture precalcolate:
% NB_Tampere_DD1C / DD5A / TUBA / C10902
% NB_Doorn_DD1C / DD5A / TUBA / C10902
% Struttura: (meaPlate).(DIV).(wellName)
% Tampere fields: starts, ends, ...
% Doorn fields: starts, ends, dur, ...
% lower = Q1 - k*IQR, upper = Q3 + k*IQR
% default k = 1.5 (conservativo: usa k=3 per outliers estremi)

%% ===== USER SETTINGS =====
selectedDIVs = {'DIV53','DIV62','DIV67','DIV83','DIV88'};
DIVs         = {'DIV53','DIV60','DIV67','DIV83','DIV90'};
saveFigures  = true;
saveFolder   = fullfile(outRoot, 'all_outputs', 'NBD_CVIBI');
outlierK     = 1.5;   % fattore IQR: 1.5 = outliers standard, 3 = solo estremi
R_refract    = 0.005; % [s] refractoriness constant per LvR (tipicamente 5 ms)

if saveFigures && ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end

%% ===== COLORI E MARKERS =====
colors = [0.70 0.10 0.10;   % DD1C  (severe)
          0.95 0.55 0.10;   % DD5A  (mild)
          0.10 0.30 0.70;   % TUBA  (ctrl)
          0.30 0.60 0.85];  % C10902 (ctrl)
markers = {'o','s','d','^'}; %#ok
clNames = {'DD1C','DD5A','TUBA','C10902'};

%% ===== LOOP SUI DUE METODI =====
methods     = {'Tampere','Doorn'};
varPrefixes = {'NB_Tampere_','NB_Doorn_'};

for mi = 1:2
    methodName = methods{mi};
    prefix     = varPrefixes{mi};

    nCL  = numel(clNames);
    nDIV = numel(selectedDIVs);

    MeanDur = nan(nCL, nDIV);
    SEM_Dur = nan(nCL, nDIV);
    N_Dur   = zeros(nCL, nDIV);
    CV_data  = cell(nCL, nDIV);
    LvR_data = cell(nCL, nDIV);   % <-- NEW

    availCL = false(1, nCL);

    %% ===== MAIN LOOP =====
    for c = 1:nCL
        varName = [prefix clNames{c}];
        if ~exist(varName,'var')
            fprintf('[%s] %s non trovato — skip\n', methodName, varName);
            continue
        end
        availCL(c) = true;
        nbData     = eval(varName);
        plateNames = fieldnames(nbData);

        for d = 1:nDIV
            DIV = selectedDIVs{d};

            dur_raw = [];
            cv_raw  = [];
            lvr_raw = [];   % <-- NEW

            for p = 1:numel(plateNames)
                meaPlate = plateNames{p};
                if ~isfield(nbData.(meaPlate), DIV), continue; end

                wellNames = fieldnames(nbData.(meaPlate).(DIV));

                for w = 1:numel(wellNames)
                    wellName = wellNames{w};
                    nb = nbData.(meaPlate).(DIV).(wellName);

                    if ~isstruct(nb) || ...
                            ~isfield(nb,'starts') || ~isfield(nb,'ends')
                        continue
                    end

                    starts = nb.starts(:);
                    ends   = nb.ends(:);
                    if isempty(starts), continue; end

                    % --- NB duration media del well ---
                    durations = ends - starts;
                    dur_raw(end+1) = mean(durations,'omitnan'); %#ok

                    % --- CV IBI e LvR del well ---
                    if numel(starts) >= 2
                        ibi = diff(sort(starts));       % inter-burst intervals
                        ibi = ibi(ibi > 0);             % rimuovi zeri

                        % CV IBI
                        mu = mean(ibi,'omitnan');
                        sg = std(ibi, 0,'omitnan');
                        if mu > 0
                            cv_raw(end+1) = sg / mu; %#ok
                        else
                            cv_raw(end+1) = NaN;     %#ok
                        end

                        % LvR --------------------------------------------------------
                        % Formula: LvR = 3/(n-1) * sum (1 - 4*Ii*Ii+1/(Ii+Ii+1)^2)
                        %                              * (1 + 4R/(Ii+Ii+1))
                        % dove n = numero di IBI, R = refractoriness constant
                        n_ibi = numel(ibi);
                        if n_ibi >= 2
                            Ii   = ibi(1:end-1);
                            Ii1  = ibi(2:end);
                            S    = Ii + Ii1;               % somma coppie
                            term = (1 - (4.*Ii.*Ii1) ./ (S.^2)) ...
                                 .* (1 + 4*R_refract ./ S);
                            lvr_val = (3 / (n_ibi - 1)) * sum(term, 'omitnan');
                            lvr_raw(end+1) = lvr_val; %#ok
                        else
                            lvr_raw(end+1) = NaN; %#ok
                        end
                    end
                end
            end

            % --- rimuovi outliers IQR ---
            dur_clean = removeOutliersIQR(dur_raw, outlierK);
            cv_clean  = removeOutliersIQR(cv_raw,  outlierK);
            lvr_clean = removeOutliersIQR(lvr_raw, outlierK);   % <-- NEW

            % --- media e 95% CI (t-student) ---
            if ~isempty(dur_clean)
                MeanDur(c,d) = mean(dur_clean,'omitnan');
                N_Dur(c,d)   = numel(dur_clean);
                n__ = numel(dur_clean);
                if n__ <= 1
                    SEM_Dur(c,d) = NaN;
                else
                    SEM_Dur(c,d) = tinv(0.975, n__-1) * std(dur_clean, 0, 'omitnan') / sqrt(n__);
                end
            end

            CV_data{c,d}  = cv_clean(:);
            LvR_data{c,d} = lvr_clean(:);   % <-- NEW
        end
    end

    idxAvail = find(availCL);
    if isempty(idxAvail)
        fprintf('[%s] Nessun dato disponibile.\n', methodName);
        continue
    end

    %% ===== FIGURA 1: NB DURATION LINE PLOT =====
    fig1 = figure('Color','w','Position',[100 100 900 560]);
    hold on; ax1 = gca;
    x = 1:nDIV;

    for c = idxAvail
        y   = MeanDur(c,:);
        sem = SEM_Dur(c,:);
        vm  = ~isnan(y) & ~isnan(sem);
        if ~any(~isnan(y)), continue; end

        % --- Shaded CI band (ci_low clampato a 0) ---
        if any(vm)
            xv = x(vm);
            yv = y(vm);
            ev = sem(vm);
            xfill = [xv, fliplr(xv)];
            yfill = [yv+ev, fliplr(max(yv-ev, 0))];
            fill(ax1, xfill, yfill, colors(c,:), 'FaceAlpha', 0.15, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        % --- Line + filled circle markers ---
        valid = ~isnan(y);
        plot(ax1, x(valid), y(valid), '-o', ...
            'Color', colors(c,:), 'LineWidth', 1.8, ...
            'MarkerSize', 6, 'MarkerFaceColor', colors(c,:), ...
            'DisplayName', clNames{c});
    end

    xlim([0.5 nDIV+0.5]);
    xticks(x); xticklabels(DIVs);
    xlabel('DIV', 'FontSize', 15);
    ylabel('Mean network burst duration (s) ± 95% CI', 'FontSize', 15);
%     title(sprintf('Network burst duration — %s (outliers removed, k=%.1f)', ...
%         methodName, outlierK), 'FontSize', 13);
    legend('Location','best','FontSize',10);
    grid on; box on;

    if saveFigures
        fname = sprintf('NB_duration_%s.png', methodName);
        exportgraphics(fig1, fullfile(saveFolder, fname), 'Resolution', 300);
        fprintf('Saved: %s\n', fname);
    end

    %% ===== FIGURA 2: CV IBI BOXPLOT =====
    fig2 = figure('Color','w','Position',[100 100 1400 600]);

    for d = 1:nDIV
        subplot(1, nDIV, d); hold on;
        ax_pos    = get(gca, 'Position');
        top_margin = 0.13;
        ax_new_h  = ax_pos(4) - top_margin + ax_pos(2) - 0.02;
        if ax_new_h > 0.1
            set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3), ax_new_h]);
        end

        dataVec  = [];
        groupVec = [];

        for ci = 1:numel(idxAvail)
            c    = idxAvail(ci);
            vals = CV_data{c,d};
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            dataVec  = [dataVec;  vals(:)];
            groupVec = [groupVec; ci*ones(numel(vals),1)];
        end

        if isempty(dataVec)
%             title(DIVs{d}, 'FontSize', 12);
            ylabel('CV IBI', 'FontSize', 15);
            xticks(1:numel(idxAvail));
            xticklabels(clNames(idxAvail));
            box on;
            continue;
        end

        clrMat = cell2mat(arrayfun(@(ci) colors(idxAvail(ci),:), ...
            1:numel(idxAvail), 'UniformOutput', false)');

        bp = boxplot(dataVec, groupVec, ...
            'Labels',  clNames(idxAvail), ...
            'Colors',  clrMat, ...
            'Symbol',  'k+', ...
            'Whisker', 1.5);
        set(bp, 'LineWidth', 1.5);

        h = findobj(gca, 'Tag', 'Box');
        for j = 1:min(numel(h), numel(idxAvail))
            patch(get(h(j),'XData'), get(h(j),'YData'), ...
                colors(idxAvail(numel(h)-j+1),:), ...
                'FaceAlpha', 0.35, ...
                'EdgeColor', colors(idxAvail(numel(h)-j+1),:), ...
                'LineWidth', 1.2);
        end

        for ci = 1:numel(idxAvail)
            c    = idxAvail(ci);
            vals = CV_data{c,d};
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            m  = mean(vals, 'omitnan');
            q3 = prctile(vals, 75);
            yL = ylim;
            dy = (yL(2) - yL(1)) * 0.04;
            if q3 + dy <= yL(2) * 0.97
                yText  = q3 + dy;
                vAlign = 'bottom';
            else
                q1    = prctile(vals, 25);
                yText  = q1 - dy;
                vAlign = 'top';
            end
            text(ci, yText, sprintf('%.2f', m), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   vAlign, ...
                'FontSize',    9, ...
                'FontWeight',  'bold', ...
                'Color',       [0.2 0.2 0.2], ...
                'BackgroundColor', 'w', ...
                'Margin',  1, ...
                'Clipping','on');
        end

%         title(DIVs{d}, 'FontSize', 12);
        ylabel('CV IBI', 'FontSize', 15);
        grid on; box on;
        yL = ylim;
        if yL(2) > 0, ylim([0, yL(2)*1.1+0.01]); end
    end

%     sgtitle(sprintf('CV Inter-Burst Interval — %s (outliers removed, k=%.1f)', ...
%         methodName, outlierK), 'FontSize', 13);

    if saveFigures
        fname = sprintf('CV_IBI_boxplot_%s.png', methodName);
        exportgraphics(fig2, fullfile(saveFolder, fname), 'Resolution', 300);
        fprintf('Saved: %s\n', fname);
    end

    %% ===== FIGURA 3: LvR BOXPLOT =====
    fig3 = figure('Color','w','Position',[100 100 1400 600]);

    for d = 1:nDIV
        subplot(1, nDIV, d); hold on;
        ax_pos    = get(gca, 'Position');
        top_margin = 0.13;
        ax_new_h  = ax_pos(4) - top_margin + ax_pos(2) - 0.02;
        if ax_new_h > 0.1
            set(gca, 'Position', [ax_pos(1), ax_pos(2), ax_pos(3), ax_new_h]);
        end

        dataVec  = [];
        groupVec = [];

        for ci = 1:numel(idxAvail)
            c    = idxAvail(ci);
            vals = LvR_data{c,d};
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            dataVec  = [dataVec;  vals(:)];
            groupVec = [groupVec; ci*ones(numel(vals),1)];
        end

        if isempty(dataVec)
%             title(DIVs{d}, 'FontSize', 12);
            ylabel('LvR', 'FontSize', 15);
            xticks(1:numel(idxAvail));
            xticklabels(clNames(idxAvail));
            box on;
            continue;
        end

        clrMat = cell2mat(arrayfun(@(ci) colors(idxAvail(ci),:), ...
            1:numel(idxAvail), 'UniformOutput', false)');

        bp = boxplot(dataVec, groupVec, ...
            'Labels',  clNames(idxAvail), ...
            'Colors',  clrMat, ...
            'Symbol',  'k+', ...
            'Whisker', 1.5);
        set(bp, 'LineWidth', 1.5);

        h = findobj(gca, 'Tag', 'Box');
        for j = 1:min(numel(h), numel(idxAvail))
            patch(get(h(j),'XData'), get(h(j),'YData'), ...
                colors(idxAvail(numel(h)-j+1),:), ...
                'FaceAlpha', 0.35, ...
                'EdgeColor', colors(idxAvail(numel(h)-j+1),:), ...
                'LineWidth', 1.2);
        end

        % --- Annotazione media sopra ogni box ---
        for ci = 1:numel(idxAvail)
            c    = idxAvail(ci);
            vals = LvR_data{c,d};
            vals = vals(~isnan(vals));
            if isempty(vals), continue; end
            m  = mean(vals, 'omitnan');
            q3 = prctile(vals, 75);
            yL = ylim;
            dy = (yL(2) - yL(1)) * 0.04;
            if q3 + dy <= yL(2) * 0.97
                yText  = q3 + dy;
                vAlign = 'bottom';
            else
                q1     = prctile(vals, 25);
                yText  = q1 - dy;
                vAlign = 'top';
            end
            text(ci, yText, sprintf('%.2f', m), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment',   vAlign, ...
                'FontSize',    9, ...
                'FontWeight',  'bold', ...
                'Color',       [0.2 0.2 0.2], ...
                'BackgroundColor', 'w', ...
                'Margin',  1, ...
                'Clipping','on');
        end

        % --- Linea di riferimento LvR = 1 (Poisson process) ---
        yL = ylim;
        plot(xlim, [1 1], '--', 'Color', [0.5 0.5 0.5], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');

%         title(DIVs{d}, 'FontSize', 12);
        ylabel(sprintf('LvR  (R = %.0f ms)', R_refract*1000), 'FontSize', 15);
        grid on; box on;
        yL = ylim;
        if yL(2) > 0, ylim([0, yL(2)*1.1+0.01]); end
    end

%     sgtitle(sprintf('LvR Inter-Burst Interval — %s (R = %.0f ms, outliers removed, k=%.1f)', ...
%         methodName, R_refract*1000, outlierK), 'FontSize', 13);

    if saveFigures
        fname = sprintf('LvR_IBI_boxplot_%s.png', methodName);
        exportgraphics(fig3, fullfile(saveFolder, fname), 'Resolution', 300);
        fprintf('Saved: %s\n', fname);
    end

    %% ===== SUMMARY CONSOLE =====
    fprintf('\n===== %s — NB DURATION, CV IBI & LvR (outliers removed, k=%.1f, R=%.0fms) =====\n', ...
        upper(methodName), outlierK, R_refract*1000);
    for c = idxAvail
        fprintf('\n  %s\n', clNames{c});
        for d = 1:nDIV
            cv_vec  = CV_data{c,d};
            cv_vec  = cv_vec(~isnan(cv_vec));
            lvr_vec = LvR_data{c,d};
            lvr_vec = lvr_vec(~isnan(lvr_vec));
            fprintf('  %s | mean_dur=%.4f s | SEM=%.4f | n_wells=%d | CV_IBI=%.3f (n=%d) | LvR=%.3f (n=%d)\n', ...
                selectedDIVs{d}, MeanDur(c,d), SEM_Dur(c,d), N_Dur(c,d), ...
                mean(cv_vec,'omitnan'), numel(cv_vec), ...
                mean(lvr_vec,'omitnan'), numel(lvr_vec));
        end
    end

end % fine loop metodi

fprintf('\nDone.\n');

end

function clean = removeOutliersIQR(data, k)
% lower fence = Q1 - k*IQR
% upper fence = Q3 + k*IQR
% k = 1.5 -> outliers standard
% k = 3.0 -> solo outliers estremi

data = data(~isnan(data));
if numel(data) < 4
    clean = data;
    return
end

q1   = prctile(data, 25);
q3   = prctile(data, 75);
iqr_ = q3 - q1;
lo   = q1 - k * iqr_;
hi   = q3 + k * iqr_;
clean = data(data >= lo & data <= hi);
end
