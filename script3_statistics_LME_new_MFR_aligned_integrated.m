function script3_statistics_LME_new_MFR_aligned_integrated()
% script3_statistics_LME_new_MFR_aligned_integrated  Full LME statistics with Plate random effect: overall tests, pairwise post-hoc (FDR), heatmaps and CI line plots.
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

% SCRIPT 3 - STATISTICAL ANALYSIS WITH PLATE AS RANDOM EFFECT (LME)
% Overall: LME con Plate come random intercept
% Post-hoc pairwise: test adattivo (normalita' -> t/Welch, altrimenti Mann-Whitney)
% Correzione multipla: FDR (Benjamini-Hochberg)
% Heatmap: -log10(p_adj) | grigio=no data, bianco=ns, blu scuro=significativo
% Line plots: mean +/- 95% CI con IQR removal (k=1.5) solo per visualizzazione
% MODIFICA: DIV 62 mostrato come "60" su tutti i grafici (DIVlabels)

%% ====================== USER PARAMETERS ================================
% Raw spikes (<CL>) and NB structs (NB_Tampere_<CL>, NB_Tampere_<CL>_V2,
% NB_Doorn_<CL>) are read directly from the BASE WORKSPACE (left there by
% run_preprocessing_processing.m) — no disk loading is needed here.
output_path  = fullfile(outRoot, 'all_outputs', 'script3');

re_output_path = fullfile(outRoot, 'all_outputs', 'script3');
fig_path       = fullfile(outRoot, 'all_outputs', 'script3', 'heatmaps');
fig_path_plots = fullfile(outRoot, 'all_outputs', 'script3', 'lineplots');

cell_lines = {'DD1C','DD5A','TUBA','C10902'};
cl_labels  = {'DD1C (severe)','DD5A (mild)','TUBA (ctrl)','C10902 (ctrl)'};
DIVs       = [53, 62, 67, 83, 88];
DIVlabels  = {'53', '60', '67', '83', '90'};  % alias visivo: DIV 62 mostrato come "60"

recording_dur_s  = 600;
recording_dur_m  = recording_dur_s/60;
active_thr_hz    = 0.1;
synch_thr_hz     = 0.2;                       % active-electrode threshold (thesis 3.3.1)
SampFreq         = 12500;                     % [Hz]
[~, ~, ParamsBurst] = get_parameter_set(2);   % single-electrode burst params (Set 2)
lvr_R            = 5e-3;                       % LvR refractoriness [s] (Shinomoto 2009; inter-NB-interval sequence, per detector)
detectors        = {'TampereV1','Doorn','TampereV2'};
alpha            = 0.05;
outlierK         = 1.5;

palette = [ 0.70 0.10 0.10;
            0.95 0.55 0.10;
            0.10 0.30 0.70;
            0.30 0.60 0.85];

det_palette = struct( ...
    'TampereV1', [0.85 0.30 0.10], ...
    'Doorn',     [0.10 0.30 0.70], ...
    'TampereV2', [0.00 0.65 0.25]);

if ~exist(output_path,'dir');    mkdir(output_path);    end
if ~exist(re_output_path,'dir'); mkdir(re_output_path); end
if ~exist(fig_path,'dir');       mkdir(fig_path);       end
if ~exist(fig_path_plots,'dir'); mkdir(fig_path_plots); end

%% ====================== LOAD DATA & COMPUTE METRICS ====================
fprintf('\n===== LOADING DATA AND COMPUTING METRICS =====\n');
T = table();

for ic = 1:numel(cell_lines)
    cl = cell_lines{ic};
    fprintf('  Processing cell line: %s\n', cl);

    % --- raw spikes from base workspace ---
    if ~evalin('base', sprintf('exist(''%s'',''var'')', cl))
        warning('Cell line %s not found in base workspace; skipping.', cl); continue;
    end
    spk = evalin('base', cl);

    NB = struct();
    for id = 1:numel(detectors)
        det = detectors{id};
        switch det
            case 'TampereV1', vname = ['NB_Tampere_' cl];
            case 'Doorn',     vname = ['NB_Doorn_'   cl];
            case 'TampereV2', vname = ['NB_Tampere_' cl '_V2'];
        end
        % --- NB structs from base workspace ---
        if ~evalin('base', sprintf('exist(''%s'',''var'')', vname))
            warning('%s not found in base workspace; skipping detector.', vname); continue;
        end
        NB.(det) = evalin('base', vname);
    end

    plates = fieldnames(spk);
    for ip = 1:numel(plates)
        plate     = plates{ip};
        div_names = fieldnames(spk.(plate));
        for idiv = 1:numel(div_names)
            div_name = div_names{idiv};
            div_num  = sscanf(div_name,'DIV%d');
            if isempty(div_num); div_num = sscanf(div_name,'%d'); end
            if ~ismember(div_num, DIVs); continue; end

            wells = fieldnames(spk.(plate).(div_name));
            for iw = 1:numel(wells)
                w   = wells{iw};
                tbl = spk.(plate).(div_name).(w);
                if ~istable(tbl) || height(tbl)==0
                    mfr = NaN;
                elseif ~ismember('Time', tbl.Properties.VariableNames) || ...
                       ~ismember('Channel', tbl.Properties.VariableNames)
                    mfr = NaN;
                else
                    % MFR aligned to mfr_burstR_NBR_CI_fixed.m
                    % total spikes in the well / recording duration (Hz)
                    mfr = compute_mfr_like_mfr_burstR_NBR(tbl, recording_dur_s);
                end

                % Activity metrics (detector-independent): SEB rate/duration,
                % % random spikes, network synchrony (active electrodes).
                [br_w, bd_w, rs_w, synch_w] = compute_activity_metrics( ...
                    tbl, SampFreq, ParamsBurst, recording_dur_s, synch_thr_hz);

                for id = 1:numel(detectors)
                    det         = detectors{id};
                    nbr         = NaN;
                    frag_per_nb = NaN;
                    nbd         = NaN;
                    cv_ibi      = NaN;
                    lvr         = NaN;
                    fragcount   = NaN;
                    fragratio   = NaN;
                    if isfield(NB,det) && isfield(NB.(det),plate) ...
                            && isfield(NB.(det).(plate),div_name) ...
                            && isfield(NB.(det).(plate).(div_name),w)
                        wellNB = NB.(det).(plate).(div_name).(w);
                        if isfield(wellNB,'starts') && ~isempty(wellNB.starts)
                            n_nb = numel(wellNB.starts);
                            nbr  = n_nb / recording_dur_m;
                            st = wellNB.starts(:);
                            if isfield(wellNB,'ends') && ~isempty(wellNB.ends)
                                nbd = mean(wellNB.ends(:) - st, 'omitnan');   % NB duration [s]
                            end
                            if numel(st) >= 2
                                ibi    = diff(sort(st));                       % inter-NB intervals [s]
                                mu_ibi = mean(ibi,'omitnan');
                                if mu_ibi > 0, cv_ibi = std(ibi,0,'omitnan')/mu_ibi; end
                                lvr = compute_lvr(ibi, lvr_R);   % Shinomoto LvR on inter-NB intervals
                            end
                            if isfield(wellNB,'num_fragments') && ~isempty(wellNB.num_fragments)
                                nf   = wellNB.num_fragments(:);
                                nf   = nf(~isnan(nf));
                                mask = nf >= 1;
                                if any(mask)
                                    frag_per_nb = sum(nf(mask)) / sum(mask);
                                end
                                if ~isempty(nf)
                                    fragcount = sum(nf);              % total fragments per well (count)
                                    fm = nf >= 1;                     % NBs with at least 1 fragment
                                    if any(fm)
                                        fragratio = sum(nf(fm)) / sum(fm);  % mean fragments/NB among NBs with >=1 fragment
                                    end
                                end
                            end
                        end
                    end
                    T = [T; table({cl},{plate},div_num,{w},{det}, ...
                        mfr,nbr,frag_per_nb,nbd,cv_ibi,lvr,fragcount,fragratio, br_w,bd_w,synch_w,rs_w, ...
                        'VariableNames', ...
                        {'CellLine','Plate','DIV','Well','Detector', ...
                         'MFR','NBR','FragPerNB','NBD','CV_IBI','LvR','FragCount','FragRatio', ...
                         'BR','BD','Synchrony','RS'})]; %#ok
                end
            end
        end
    end
end

if isempty(T), error('No data loaded.'); end
writetable(T, fullfile(re_output_path,'metrics_per_well.xlsx'));
fprintf('  Saved metrics_per_well.xlsx (%d rows)\n', height(T));

%% ====================== ANALYSIS 1 (LME) ===============================
fprintf('\n===== ANALYSIS 1 (LME): CELL LINE x DIV (Tampere V1) =====\n');

T1               = T(strcmp(T.Detector,'TampereV1'), :);
metrics_A1       = {'MFR','NBR','FragPerNB'};
metric_labels_A1 = {'MFR (Hz)','NBR (bursts/min)','Fragments / NB'};
A1_overall       = table();
A1_pairwise      = table();

for im = 1:numel(metrics_A1)
    mname = metrics_A1{im};
    fprintf('  Metric: %s\n', mname);
    sub = T1(~isnan(T1.(mname)), :);

    [overall_pval, test_name, eff_size, posthoc_lme] = run_overall_lme( ...
        sub, mname, {'CellLine','DIV'}, 'Plate', cell_lines, DIVs);

    A1_overall = [A1_overall; table({mname},{test_name},overall_pval,eff_size, ...
        'VariableNames',{'Metric','Test','p_value','EffectSize'})]; %#ok

    pairs = nchoosek(1:numel(cell_lines),2);
    for idiv = 1:numel(DIVs)
        d = DIVs(idiv);
        for ip = 1:size(pairs,1)
            A  = cell_lines{pairs(ip,1)};
            B  = cell_lines{pairs(ip,2)};
            xA = sub.(mname)(strcmp(sub.CellLine,A) & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.CellLine,B) & sub.DIV==d);
            [p, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A1_pairwise = [A1_pairwise; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                p,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                  'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
        end
    end

    if ~isempty(posthoc_lme)
        writetable(posthoc_lme, fullfile(re_output_path, ...
            sprintf('Analysis1_LME_posthoc_%s.xlsx', mname)));
    end
end

A1_pairwise.p_adj = nan(height(A1_pairwise),1);
for im = 1:numel(metrics_A1)
    mask = strcmp(A1_pairwise.Metric, metrics_A1{im});
    A1_pairwise.p_adj(mask) = fdr_bh(A1_pairwise.p_raw(mask));
end
A1_pairwise.Significance = sig_stars(A1_pairwise.p_adj);

writetable(A1_overall,  fullfile(re_output_path,'Analysis1_overall.xlsx'));
writetable(A1_pairwise, fullfile(re_output_path,'Analysis1_pairwise.xlsx'));

% --- Heatmaps A1 ---
for im = 1:numel(metrics_A1)
    mname = metrics_A1{im};
    plot_pairwise_heatmap(A1_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 1 (LME) - %s (Tampere V1)', mname), ...
        fullfile(fig_path, sprintf('A1_heatmap_%s', mname)));
end

% --- Line plots A1 ---
fprintf('\n===== ANALYSIS 1: LINE PLOTS =====\n');
overall_pairs_A1 = compute_overall_cellline_pairs(T1, metrics_A1, cell_lines);

for im = 1:numel(metrics_A1)
    mname = metrics_A1{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T1, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_, '-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A1{im});
%     title(ax, sprintf('Analysis 1 - %s (Tampere V1)', mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A1.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A1_lineplot_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

%% ============ ANALYSIS 1B (LME): ACTIVITY METRICS (CellLine x DIV) ======
fprintf('\n===== ANALYSIS 1B (LME): ACTIVITY METRICS (CellLine x DIV) =====\n');

metrics_A1B       = {'MFR','BR','BD','Synchrony','RS'};
metric_labels_A1B = {'MFR (Hz)','SEB rate (bursts/min)','SEB duration (s)', ...
                     'Active electrodes','% Random spikes'};
A1B_overall  = table();
A1B_pairwise = table();

for im = 1:numel(metrics_A1B)
    mname = metrics_A1B{im};
    fprintf('  Metric: %s\n', mname);
    sub = T1(~isnan(T1.(mname)), :);

    [overall_pval, test_name, eff_size, posthoc_lme] = run_overall_lme( ...
        sub, mname, {'CellLine','DIV'}, 'Plate', cell_lines, DIVs);

    A1B_overall = [A1B_overall; table({mname},{test_name},overall_pval,eff_size, ...
        'VariableNames',{'Metric','Test','p_value','EffectSize'})]; %#ok

    pairs = nchoosek(1:numel(cell_lines),2);
    for idiv = 1:numel(DIVs)
        d = DIVs(idiv);
        for ip = 1:size(pairs,1)
            A  = cell_lines{pairs(ip,1)};
            B  = cell_lines{pairs(ip,2)};
            xA = sub.(mname)(strcmp(sub.CellLine,A) & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.CellLine,B) & sub.DIV==d);
            [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A1B_pairwise = [A1B_pairwise; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                  'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
        end
    end

    if ~isempty(posthoc_lme)
        writetable(posthoc_lme, fullfile(re_output_path, ...
            sprintf('Analysis1B_LME_posthoc_%s.xlsx', mname)));
    end
end

A1B_pairwise.p_adj = nan(height(A1B_pairwise),1);
for im = 1:numel(metrics_A1B)
    mask = strcmp(A1B_pairwise.Metric, metrics_A1B{im});
    A1B_pairwise.p_adj(mask) = fdr_bh(A1B_pairwise.p_raw(mask));
end
A1B_pairwise.Significance = sig_stars(A1B_pairwise.p_adj);

writetable(A1B_overall,  fullfile(re_output_path,'Analysis1B_activity_overall.xlsx'));
writetable(A1B_pairwise, fullfile(re_output_path,'Analysis1B_activity_pairwise.xlsx'));

% --- Heatmaps A1B ---
for im = 1:numel(metrics_A1B)
    mname = metrics_A1B{im};
    plot_pairwise_heatmap(A1B_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 1B (LME) - %s', mname), ...
        fullfile(fig_path, sprintf('A1B_heatmap_%s', mname)));
end

% --- Line plots A1B ---
fprintf('\n===== ANALYSIS 1B: LINE PLOTS =====\n');
overall_pairs_A1B = compute_overall_cellline_pairs(T1, metrics_A1B, cell_lines);

for im = 1:numel(metrics_A1B)
    mname = metrics_A1B{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T1, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_, '-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A1B{im});
%     title(ax, sprintf('Analysis 1B - %s', mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A1B.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A1B_lineplot_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

%% ============ ANALYSIS 1C (LME): NB METRICS (CellLine x DIV, Tampere V1) =
fprintf('\n===== ANALYSIS 1C (LME): NB METRICS (CellLine x DIV, Tampere V1) =====\n');

metrics_A1C       = {'NBR','NBD','CV_IBI','LvR'};
metric_labels_A1C = {'NBR (bursts/min)','NB duration (s)','CV(IBI)','LvR'};
A1C_overall  = table();
A1C_pairwise = table();

for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    fprintf('  Metric: %s\n', mname);
    sub = T1(~isnan(T1.(mname)), :);

    [overall_pval, test_name, eff_size, posthoc_lme] = run_overall_lme( ...
        sub, mname, {'CellLine','DIV'}, 'Plate', cell_lines, DIVs);

    A1C_overall = [A1C_overall; table({mname},{test_name},overall_pval,eff_size, ...
        'VariableNames',{'Metric','Test','p_value','EffectSize'})]; %#ok

    pairs = nchoosek(1:numel(cell_lines),2);
    for idiv = 1:numel(DIVs)
        d = DIVs(idiv);
        for ip = 1:size(pairs,1)
            A  = cell_lines{pairs(ip,1)};
            B  = cell_lines{pairs(ip,2)};
            xA = sub.(mname)(strcmp(sub.CellLine,A) & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.CellLine,B) & sub.DIV==d);
            [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A1C_pairwise = [A1C_pairwise; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                  'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
        end
    end

    if ~isempty(posthoc_lme)
        writetable(posthoc_lme, fullfile(re_output_path, ...
            sprintf('Analysis1C_LME_posthoc_%s.xlsx', mname)));
    end
end

A1C_pairwise.p_adj = nan(height(A1C_pairwise),1);
for im = 1:numel(metrics_A1C)
    mask = strcmp(A1C_pairwise.Metric, metrics_A1C{im});
    A1C_pairwise.p_adj(mask) = fdr_bh(A1C_pairwise.p_raw(mask));
end
A1C_pairwise.Significance = sig_stars(A1C_pairwise.p_adj);

writetable(A1C_overall,  fullfile(re_output_path,'Analysis1C_NBmetrics_overall.xlsx'));
writetable(A1C_pairwise, fullfile(re_output_path,'Analysis1C_NBmetrics_pairwise.xlsx'));

% --- Heatmaps A1C ---
for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    plot_pairwise_heatmap(A1C_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 1C (LME) - %s (Tampere V1)', mname), ...
        fullfile(fig_path, sprintf('A1C_heatmap_%s', mname)));
end

% --- Line plots A1C ---
fprintf('\n===== ANALYSIS 1C: LINE PLOTS =====\n');
overall_pairs_A1C = compute_overall_cellline_pairs(T1, metrics_A1C, cell_lines);

for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T1, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_, '-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A1C{im});
%     title(ax, sprintf('Analysis 1C - %s (Tampere V1)', mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A1C.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A1C_lineplot_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

%% ====== ANALYSIS 1C-Doorn (LME): NB METRICS (CellLine x DIV, Doorn) ======
% Same cross-cell-line LME as Analysis 1C (NBR, NBD, CV_IBI, LvR), but using
% the Doorn detector instead of Tampere V1. Outputs use an "A1C_Doorn"/
% "Analysis1C_Doorn" prefix so nothing overwrites the Tampere V1 results.
fprintf('\n===== ANALYSIS 1C-Doorn (LME): NB METRICS (CellLine x DIV, Doorn) =====\n');

T1_doorn      = T(strcmp(T.Detector,'Doorn'), :);
A1Cd_overall  = table();
A1Cd_pairwise = table();

for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    fprintf('  Metric: %s\n', mname);
    sub = T1_doorn(~isnan(T1_doorn.(mname)), :);

    [overall_pval, test_name, eff_size, posthoc_lme] = run_overall_lme( ...
        sub, mname, {'CellLine','DIV'}, 'Plate', cell_lines, DIVs);

    A1Cd_overall = [A1Cd_overall; table({mname},{test_name},overall_pval,eff_size, ...
        'VariableNames',{'Metric','Test','p_value','EffectSize'})]; %#ok

    pairs = nchoosek(1:numel(cell_lines),2);
    for idiv = 1:numel(DIVs)
        d = DIVs(idiv);
        for ip = 1:size(pairs,1)
            A  = cell_lines{pairs(ip,1)};
            B  = cell_lines{pairs(ip,2)};
            xA = sub.(mname)(strcmp(sub.CellLine,A) & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.CellLine,B) & sub.DIV==d);
            [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A1Cd_pairwise = [A1Cd_pairwise; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                  'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
        end
    end

    if ~isempty(posthoc_lme)
        writetable(posthoc_lme, fullfile(re_output_path, ...
            sprintf('Analysis1C_Doorn_LME_posthoc_%s.xlsx', mname)));
    end
end

A1Cd_pairwise.p_adj = nan(height(A1Cd_pairwise),1);
for im = 1:numel(metrics_A1C)
    mask = strcmp(A1Cd_pairwise.Metric, metrics_A1C{im});
    A1Cd_pairwise.p_adj(mask) = fdr_bh(A1Cd_pairwise.p_raw(mask));
end
A1Cd_pairwise.Significance = sig_stars(A1Cd_pairwise.p_adj);

writetable(A1Cd_overall,  fullfile(re_output_path,'Analysis1C_Doorn_NBmetrics_overall.xlsx'));
writetable(A1Cd_pairwise, fullfile(re_output_path,'Analysis1C_Doorn_NBmetrics_pairwise.xlsx'));

% --- Heatmaps A1C-Doorn ---
for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    plot_pairwise_heatmap(A1Cd_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 1C (LME) - %s (Doorn)', mname), ...
        fullfile(fig_path, sprintf('A1C_Doorn_heatmap_%s', mname)));
end

% --- Line plots A1C-Doorn ---
fprintf('\n===== ANALYSIS 1C-Doorn: LINE PLOTS =====\n');
overall_pairs_A1Cd = compute_overall_cellline_pairs(T1_doorn, metrics_A1C, cell_lines);

for im = 1:numel(metrics_A1C)
    mname = metrics_A1C{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T1_doorn, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_, '-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A1C{im});
%     title(ax, sprintf('Analysis 1C - %s (Doorn)', mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A1Cd.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A1C_Doorn_lineplot_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

%% ====== ANALYSIS 1D (LME): FRAGMENTATION (CellLine x DIV) ================
% Cross-cell-line LME (same design as Analysis 1C) for the fragmentation
% metrics FragPerNB and FragCount, computed separately for the Tampere V1
% and Doorn detectors. Outputs are tagged "A1D_<detector>" so nothing
% collides with the other analyses.
fprintf('\n===== ANALYSIS 1D (LME): FRAGMENTATION (CellLine x DIV) =====\n');

metrics_A1D       = {'FragPerNB','FragCount'};
metric_labels_A1D = {'Fragments / NB','Fragment count'};
frag1D_dets       = {'TampereV1','Doorn'};

for jdet = 1:numel(frag1D_dets)
    detTag = frag1D_dets{jdet};
    Td     = T(strcmp(T.Detector, detTag), :);

    A1D_overall  = table();
    A1D_pairwise = table();

    for im = 1:numel(metrics_A1D)
        mname = metrics_A1D{im};
        fprintf('  [%s] Metric: %s\n', detTag, mname);
        sub = Td(~isnan(Td.(mname)), :);

        [overall_pval, test_name, eff_size, posthoc_lme] = run_overall_lme( ...
            sub, mname, {'CellLine','DIV'}, 'Plate', cell_lines, DIVs);

        A1D_overall = [A1D_overall; table({mname},{test_name},overall_pval,eff_size, ...
            'VariableNames',{'Metric','Test','p_value','EffectSize'})]; %#ok

        pairs = nchoosek(1:numel(cell_lines),2);
        for idiv = 1:numel(DIVs)
            d = DIVs(idiv);
            for ip = 1:size(pairs,1)
                A  = cell_lines{pairs(ip,1)};
                B  = cell_lines{pairs(ip,2)};
                xA = sub.(mname)(strcmp(sub.CellLine,A) & sub.DIV==d);
                xB = sub.(mname)(strcmp(sub.CellLine,B) & sub.DIV==d);
                [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
                A1D_pairwise = [A1D_pairwise; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                    pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                    'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                      'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
            end
        end

        if ~isempty(posthoc_lme)
            writetable(posthoc_lme, fullfile(re_output_path, ...
                sprintf('Analysis1D_%s_LME_posthoc_%s.xlsx', detTag, mname)));
        end
    end

    A1D_pairwise.p_adj = nan(height(A1D_pairwise),1);
    for im = 1:numel(metrics_A1D)
        mask = strcmp(A1D_pairwise.Metric, metrics_A1D{im});
        A1D_pairwise.p_adj(mask) = fdr_bh(A1D_pairwise.p_raw(mask));
    end
    A1D_pairwise.Significance = sig_stars(A1D_pairwise.p_adj);

    writetable(A1D_overall,  fullfile(re_output_path, sprintf('Analysis1D_%s_Frag_overall.xlsx',  detTag)));
    writetable(A1D_pairwise, fullfile(re_output_path, sprintf('Analysis1D_%s_Frag_pairwise.xlsx', detTag)));

    % --- Heatmaps A1D ---
    for im = 1:numel(metrics_A1D)
        mname = metrics_A1D{im};
        plot_pairwise_heatmap(A1D_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
            sprintf('Analysis 1D (LME) - %s (%s)', mname, detTag), ...
            fullfile(fig_path, sprintf('A1D_%s_heatmap_%s', detTag, mname)));
    end

    % --- Line plots A1D ---
    overall_pairs_A1D = compute_overall_cellline_pairs(Td, metrics_A1D, cell_lines);
    for im = 1:numel(metrics_A1D)
        mname = metrics_A1D{im};
        fig = figure('Color','w','Position',[100 100 950 500]);
        ax  = axes('Position',[0.08 0.12 0.65 0.78]);
        hold(ax,'on');

        handles = gobjects(numel(cell_lines),1);
        for ic = 1:numel(cell_lines)
            cl = cell_lines{ic};
            [means_, ci_low, ci_up] = mean_ci_per_div(Td, cl, mname, DIVs, alpha, outlierK);
            valid = ~isnan(means_);
            if any(valid)
                xfill = [DIVs(valid), fliplr(DIVs(valid))];
                yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
                fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                    'EdgeColor','none','HandleVisibility','off');
            end
            handles(ic) = plot(ax, DIVs, means_, '-o', ...
                'Color',palette(ic,:), 'LineWidth',1.8, ...
                'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
        end

        xlabel(ax,'DIV'); ylabel(ax, metric_labels_A1D{im});
%         title(ax, sprintf('Analysis 1D - %s (%s)', mname, detTag),'FontWeight','bold');
        place_legend_outside(legend(ax, handles, cl_labels));
        set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
        grid(ax,'on'); ax.GridAlpha = 0.15;

        draw_overall_brackets_cellline(fig, ax, overall_pairs_A1D.(mname), cell_lines, palette, DIVs);

        savebase = fullfile(fig_path_plots, sprintf('A1D_%s_lineplot_%s', detTag, mname));
        exportgraphics(fig, [savebase '.png'],'Resolution',300);
        exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
        close(fig);
    end
end

%% ====================== ANALYSIS 2 (LME) ===============================
fprintf('\n===== ANALYSIS 2 (LME): Tampere V1 vs Doorn =====\n');

metrics_A23       = {'NBR','FragPerNB'};
metric_labels_A23 = {'NBR (bursts/min)','Fragments / NB'};
A2_overall        = table();
A2_pairwise       = table();

for ic = 1:numel(cell_lines)
    cl = cell_lines{ic};
    for im = 1:numel(metrics_A23)
        mname = metrics_A23{im};
        sub = T(strcmp(T.CellLine,cl) & ismember(T.Detector,{'TampereV1','Doorn'}) ...
              & ~isnan(T.(mname)), :);

        [p_ov, tname_ov, d_ov] = run_overall_lme_detector(sub, mname, 'Plate', 'TampereV1','Doorn');

        A2_overall = [A2_overall; table({cl},{mname}, ...
            sum(strcmp(sub.Detector,'TampereV1')), ...
            sum(strcmp(sub.Detector,'Doorn')), ...
            p_ov,{tname_ov},d_ov,'VariableNames', ...
            {'CellLine','Metric','n_V1','n_Doorn','p_raw','Test','Cohen_d'})]; %#ok

        for idiv = 1:numel(DIVs)
            d  = DIVs(idiv);
            xA = sub.(mname)(strcmp(sub.Detector,'TampereV1') & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.Detector,'Doorn')     & sub.DIV==d);
            [p, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A2_pairwise = [A2_pairwise; table({cl},{mname},d,numel(xA),numel(xB), ...
                p,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames', ...
                {'CellLine','Metric','DIV','n_V1','n_Doorn','p_raw','Test','Cohen_d','MeanDiff_V1minusDoorn'})]; %#ok
        end
    end
end

A2_overall.p_adj        = fdr_bh(A2_overall.p_raw);
A2_overall.Significance = sig_stars(A2_overall.p_adj);
A2_pairwise.p_adj       = nan(height(A2_pairwise),1);
for im = 1:numel(metrics_A23)
    mask = strcmp(A2_pairwise.Metric, metrics_A23{im});
    A2_pairwise.p_adj(mask) = fdr_bh(A2_pairwise.p_raw(mask));
end
A2_pairwise.Significance = sig_stars(A2_pairwise.p_adj);

writetable(A2_overall,  fullfile(re_output_path,'Analysis2_overall.xlsx'));
writetable(A2_pairwise, fullfile(re_output_path,'Analysis2_pairwise.xlsx'));

% --- Heatmaps A2 ---
for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    plot_detector_heatmap(A2_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 2 (LME) - %s (V1 vs Doorn)', mname), ...
        fullfile(fig_path, sprintf('A2_heatmap_%s', mname)), 'V1','Doorn');
end

% --- Analysis 2: extra heatmaps (V1 vs Doorn) for NB-timing metrics ---
metrics_A2H = {'NBR','NBD','CV_IBI','LvR'};
A2H_pairwise = table();
for ic = 1:numel(cell_lines)
    cl = cell_lines{ic};
    for im = 1:numel(metrics_A2H)
        mname = metrics_A2H{im};
        sub = T(strcmp(T.CellLine,cl) & ismember(T.Detector,{'TampereV1','Doorn'}) ...
              & ~isnan(T.(mname)), :);
        for idiv = 1:numel(DIVs)
            d  = DIVs(idiv);
            xA = sub.(mname)(strcmp(sub.Detector,'TampereV1') & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.Detector,'Doorn')     & sub.DIV==d);
            [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A2H_pairwise = [A2H_pairwise; table({cl},{mname},d,numel(xA),numel(xB), ...
                pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames', ...
                {'CellLine','Metric','DIV','n_V1','n_Doorn','p_raw','Test','Cohen_d','MeanDiff_V1minusDoorn'})]; %#ok
        end
    end
end
A2H_pairwise.p_adj = nan(height(A2H_pairwise),1);
for im = 1:numel(metrics_A2H)
    mask = strcmp(A2H_pairwise.Metric, metrics_A2H{im});
    A2H_pairwise.p_adj(mask) = fdr_bh(A2H_pairwise.p_raw(mask));
end
A2H_pairwise.Significance = sig_stars(A2H_pairwise.p_adj);
writetable(A2H_pairwise, fullfile(re_output_path,'Analysis2_NBmetrics_pairwise.xlsx'));
for im = 1:numel(metrics_A2H)
    mname = metrics_A2H{im};
    plot_detector_heatmap(A2H_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 2 (LME) - %s (V1 vs Doorn)', mname), ...
        fullfile(fig_path, sprintf('A2_heatmap_%s', mname)), 'V1','Doorn');
end

% --- Line plots A2 - Group A: per cell line, V1 vs Doorn ---
fprintf('\n===== ANALYSIS 2: LINE PLOTS =====\n');
for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    for ic = 1:numel(cell_lines)
        cl  = cell_lines{ic};
        fig = figure('Color','w','Position',[100 100 850 450]);
        ax  = axes('Position',[0.10 0.12 0.65 0.78]);
        hold(ax,'on');

        h_v1 = plot_detector_trace(ax, T, cl, 'TampereV1', mname, DIVs, det_palette.TampereV1, alpha, outlierK);
        h_do = plot_detector_trace(ax, T, cl, 'Doorn',     mname, DIVs, det_palette.Doorn,     alpha, outlierK);

        xlabel(ax,'DIV'); ylabel(ax, metric_labels_A23{im});
%         title(ax, sprintf('Analysis 2 - %s (%s): V1 vs Doorn', cl, mname),'FontWeight','bold');
        place_legend_outside(legend(ax,[h_v1 h_do],{'Tampere V1','Doorn'}));
        set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
        grid(ax,'on'); ax.GridAlpha = 0.15;

        row = strcmp(A2_overall.CellLine,cl) & strcmp(A2_overall.Metric,mname);
        if any(row) && A2_overall.p_adj(row) < alpha
            draw_overall_bracket_detector(fig, ax, A2_overall.p_adj(row), ...
                'Tampere V1','Doorn', det_palette.TampereV1, det_palette.Doorn);
        end

        savebase = fullfile(fig_path_plots, sprintf('A2_lineplot_%s_%s', cl, mname));
        exportgraphics(fig, [savebase '.png'],'Resolution',300);
        exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
        close(fig);
    end
end

% --- Line plots A2 - Group B: 4 cell lines Doorn only ---
T_do = T(strcmp(T.Detector,'Doorn'), :);
overall_pairs_A2_B = compute_overall_cellline_pairs(T_do, metrics_A23, cell_lines);

for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T_do, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:),'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_,'-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A23{im});
%     title(ax, sprintf('Analysis 2 - %s (Doorn): cell lines comparison',mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A2_B.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A2_lineplot_Doorn_allCL_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

%% ====================== ANALYSIS 3 (LME) ===============================
fprintf('\n===== ANALYSIS 3 (LME): Tampere V1 vs Tampere V2 =====\n');

A3_overall  = table();
A3_pairwise = table();

for ic = 1:numel(cell_lines)
    cl = cell_lines{ic};
    for im = 1:numel(metrics_A23)
        mname = metrics_A23{im};
        sub = T(strcmp(T.CellLine,cl) & ismember(T.Detector,{'TampereV1','TampereV2'}) ...
              & ~isnan(T.(mname)), :);

        [p_ov, tname_ov, d_ov] = run_overall_lme_detector(sub, mname, 'Plate', 'TampereV1','TampereV2');

        A3_overall = [A3_overall; table({cl},{mname}, ...
            sum(strcmp(sub.Detector,'TampereV1')), ...
            sum(strcmp(sub.Detector,'TampereV2')), ...
            p_ov,{tname_ov},d_ov,'VariableNames', ...
            {'CellLine','Metric','n_V1','n_V2','p_raw','Test','Cohen_d'})]; %#ok

        for idiv = 1:numel(DIVs)
            d  = DIVs(idiv);
            xA = sub.(mname)(strcmp(sub.Detector,'TampereV1') & sub.DIV==d);
            xB = sub.(mname)(strcmp(sub.Detector,'TampereV2') & sub.DIV==d);
            [p, tname, d_eff] = run_pairwise_test(xA, xB, false);
            A3_pairwise = [A3_pairwise; table({cl},{mname},d,numel(xA),numel(xB), ...
                p,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                'VariableNames', ...
                {'CellLine','Metric','DIV','n_V1','n_V2','p_raw','Test','Cohen_d','MeanDiff_V1minusV2'})]; %#ok
        end
    end
end

A3_overall.p_adj        = fdr_bh(A3_overall.p_raw);
A3_overall.Significance = sig_stars(A3_overall.p_adj);
A3_pairwise.p_adj       = nan(height(A3_pairwise),1);
for im = 1:numel(metrics_A23)
    mask = strcmp(A3_pairwise.Metric, metrics_A23{im});
    A3_pairwise.p_adj(mask) = fdr_bh(A3_pairwise.p_raw(mask));
end
A3_pairwise.Significance = sig_stars(A3_pairwise.p_adj);

writetable(A3_overall,  fullfile(re_output_path,'Analysis3_overall.xlsx'));
writetable(A3_pairwise, fullfile(re_output_path,'Analysis3_pairwise.xlsx'));

% --- Heatmaps A3 ---
for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    plot_detector_heatmap(A3_pairwise, mname, cell_lines, DIVs, DIVlabels, ...
        sprintf('Analysis 3 (LME) - %s (V1 vs V2)', mname), ...
        fullfile(fig_path, sprintf('A3_heatmap_%s', mname)), 'V1','V2');
end

% --- Line plots A3 - Group A: per cell line, V1 vs V2 ---
fprintf('\n===== ANALYSIS 3: LINE PLOTS =====\n');
for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    for ic = 1:numel(cell_lines)
        cl  = cell_lines{ic};
        fig = figure('Color','w','Position',[100 100 850 450]);
        ax  = axes('Position',[0.10 0.12 0.65 0.78]);
        hold(ax,'on');

        h_v1 = plot_detector_trace(ax, T, cl, 'TampereV1', mname, DIVs, det_palette.TampereV1, alpha, outlierK);
        h_v2 = plot_detector_trace(ax, T, cl, 'TampereV2', mname, DIVs, det_palette.TampereV2, alpha, outlierK);

        xlabel(ax,'DIV'); ylabel(ax, metric_labels_A23{im});
%         title(ax, sprintf('Analysis 3 - %s (%s): V1 vs V2', cl, mname),'FontWeight','bold');
        place_legend_outside(legend(ax,[h_v1 h_v2],{'Tampere V1','Tampere V2'}));
        set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
        grid(ax,'on'); ax.GridAlpha = 0.15;

        row = strcmp(A3_overall.CellLine,cl) & strcmp(A3_overall.Metric,mname);
        if any(row) && A3_overall.p_adj(row) < alpha
            draw_overall_bracket_detector(fig, ax, A3_overall.p_adj(row), ...
                'Tampere V1','Tampere V2', det_palette.TampereV1, det_palette.TampereV2);
        end

        savebase = fullfile(fig_path_plots, sprintf('A3_lineplot_%s_%s', cl, mname));
        exportgraphics(fig, [savebase '.png'],'Resolution',300);
        exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
        close(fig);
    end
end

% --- Line plots A3 - Group B: 4 cell lines TampereV2 only ---
T_v2 = T(strcmp(T.Detector,'TampereV2'), :);
overall_pairs_A3_B = compute_overall_cellline_pairs(T_v2, metrics_A23, cell_lines);

for im = 1:numel(metrics_A23)
    mname = metrics_A23{im};
    fig = figure('Color','w','Position',[100 100 950 500]);
    ax  = axes('Position',[0.08 0.12 0.65 0.78]);
    hold(ax,'on');

    handles = gobjects(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        cl = cell_lines{ic};
        [means_, ci_low, ci_up] = mean_ci_per_div(T_v2, cl, mname, DIVs, alpha, outlierK);
        valid = ~isnan(means_);
        if any(valid)
            xfill = [DIVs(valid), fliplr(DIVs(valid))];
            yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
            fill(ax, xfill, yfill, palette(ic,:),'FaceAlpha',0.15, ...
                'EdgeColor','none','HandleVisibility','off');
        end
        handles(ic) = plot(ax, DIVs, means_,'-o', ...
            'Color',palette(ic,:), 'LineWidth',1.8, ...
            'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',palette(ic,:));
    end

    xlabel(ax,'DIV'); ylabel(ax, metric_labels_A23{im});
%     title(ax, sprintf('Analysis 3 - %s (Tampere V2): cell lines comparison',mname),'FontWeight','bold');
    place_legend_outside(legend(ax, handles, cl_labels));
    set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
    grid(ax,'on'); ax.GridAlpha = 0.15;

    draw_overall_brackets_cellline(fig, ax, overall_pairs_A3_B.(mname), cell_lines, palette, DIVs);

    savebase = fullfile(fig_path_plots, sprintf('A3_lineplot_V2_allCL_%s', mname));
    exportgraphics(fig, [savebase '.png'],'Resolution',300);
    exportgraphics(fig, [savebase '.pdf'],'ContentType','vector');
    close(fig);
end

save(fullfile(re_output_path,'all_results.mat'), 'T', ...
    'A1_overall','A1_pairwise','A2_overall','A2_pairwise','A3_overall','A3_pairwise');

%% ============ ANALYSIS 3: V1 vs V2 GRID (NBR & Frag per NB) =============
fprintf('\n===== ANALYSIS 3: V1 vs V2 SUMMARY GRID (NBR & Frag per NB) =====\n');

v3_metrics = {'NBR','FragPerNB'};
v3_labels  = {'NB rate (bursts/min)','Fragments / NB'};
colV1 = det_palette.TampereV1;
colV2 = det_palette.TampereV2;

% Overall V1-vs-V2 p per (cell line, metric), FDR within metric
pV3 = nan(numel(cell_lines), numel(v3_metrics));
for im = 1:numel(v3_metrics)
    praw = nan(numel(cell_lines),1);
    for ic = 1:numel(cell_lines)
        sub = T(strcmp(T.CellLine,cell_lines{ic}) & ...
                ismember(T.Detector,{'TampereV1','TampereV2'}) & ~isnan(T.(v3_metrics{im})), :);
        [p_ov,~,~] = run_overall_lme_detector(sub, v3_metrics{im}, 'Plate','TampereV1','TampereV2');
        praw(ic) = p_ov;
    end
    pV3(:,im) = fdr_bh(praw);
end

figV = figure('Color','w','Name','A3_V1vsV2_grid','Position',[60 40 1000 1250]);
tlV  = tiledlayout(figV, numel(cell_lines), numel(v3_metrics), ...
    'TileSpacing','compact','Padding','compact');
hV1 = []; hV2 = [];
for ic = 1:numel(cell_lines)
    cl = cell_lines{ic};
    for im = 1:numel(v3_metrics)
        mname = v3_metrics{im};
        ax = nexttile(tlV); hold(ax,'on');
        h1 = plot_detector_trace(ax, T, cl, 'TampereV1', mname, DIVs, colV1, alpha, outlierK);
        h2 = plot_detector_trace(ax, T, cl, 'TampereV2', mname, DIVs, colV2, alpha, outlierK);
        if ic==1 && im==1, hV1 = h1; hV2 = h2; end
        if ~isnan(pV3(ic,im)) && pV3(ic,im) < 0.05
            draw_overall_bracket_detector(figV, ax, pV3(ic,im), 'V1','V2', colV1, colV2);
        end
        set(ax,'FontSize',14,'FontWeight','bold','Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
        grid(ax,'on'); ax.GridAlpha = 0.15;
        if im==1
            ylabel(ax, sprintf('%s | %s', cl, v3_labels{im}), 'FontSize',15,'FontWeight','bold');
        else
            ylabel(ax, v3_labels{im}, 'FontSize',15,'FontWeight','bold');
        end
        if ic==numel(cell_lines), xlabel(ax,'DIV','FontSize',15,'FontWeight','bold'); end
    end
end
if ~isempty(hV1)
    lgd = legend([hV1 hV2], {'Tampere V1','Tampere V2'}, 'Box','off','FontSize',10);
    lgd.Layout.Tile = 'east';
end
savebaseV = fullfile(fig_path_plots,'A3_V1vsV2_grid_NBR_FragPerNB');
exportgraphics(figV, [savebaseV '.png'],'Resolution',300);
exportgraphics(figV, [savebaseV '.pdf'],'ContentType','vector');
close(figV);

%% ============ ANALYSIS 4: FRAGMENTATION (CellLine x DIV, per detector) ===
fprintf('\n===== ANALYSIS 4: FRAGMENTATION (Fragments count & ratio) =====\n');

frag_dets    = {'TampereV1','Doorn'};
frag_dtitle  = {'Tampere','Doorn'};
frag_metrics = {'FragCount','FragRatio'};
frag_mlabel  = {'Fragments count','Fragments / NB'};

A4_pairwise = table();
ovp_det     = cell(1,numel(frag_dets));

for jd = 1:numel(frag_dets)
    det  = frag_dets{jd};
    Tdet = T(strcmp(T.Detector,det), :);

    % Overall cell-line pairwise (pooled) -> significance brackets
    ovp_det{jd} = compute_overall_cellline_pairs(Tdet, frag_metrics, cell_lines);

    % Per-DIV cell-line pairwise -> heatmaps
    Pdet  = table();
    pairs = nchoosek(1:numel(cell_lines),2);
    for im = 1:numel(frag_metrics)
        mname = frag_metrics{im};
        for idiv = 1:numel(DIVs)
            d = DIVs(idiv);
            for ip = 1:size(pairs,1)
                A  = cell_lines{pairs(ip,1)};
                B  = cell_lines{pairs(ip,2)};
                xA = Tdet.(mname)(strcmp(Tdet.CellLine,A) & Tdet.DIV==d);
                xB = Tdet.(mname)(strcmp(Tdet.CellLine,B) & Tdet.DIV==d);
                [pv, tname, d_eff] = run_pairwise_test(xA, xB, false);
                Pdet = [Pdet; table({mname},d,{A},{B},numel(xA),numel(xB), ...
                    pv,{tname},d_eff,mean(xA,'omitnan')-mean(xB,'omitnan'), ...
                    'VariableNames',{'Metric','DIV','GroupA','GroupB','nA','nB', ...
                                      'p_raw','Test','Cohen_d','MeanDiff_AminusB'})]; %#ok
            end
        end
    end
    Pdet.p_adj = nan(height(Pdet),1);
    for im = 1:numel(frag_metrics)
        mask = strcmp(Pdet.Metric, frag_metrics{im});
        Pdet.p_adj(mask) = fdr_bh(Pdet.p_raw(mask));
    end
    Pdet.Significance = sig_stars(Pdet.p_adj);
    Pdet.Detector = repmat(string(det), height(Pdet), 1);
    A4_pairwise = [A4_pairwise; Pdet]; %#ok

    % Heatmaps: one per metric for this detector
    for im = 1:numel(frag_metrics)
        mname = frag_metrics{im};
        plot_pairwise_heatmap(Pdet, mname, cell_lines, DIVs, DIVlabels, ...
            sprintf('Analysis 4 - %s (%s)', frag_mlabel{im}, frag_dtitle{jd}), ...
            fullfile(fig_path, sprintf('A4_heatmap_%s_%s', mname, frag_dtitle{jd})));
    end
end
writetable(A4_pairwise, fullfile(re_output_path,'Analysis4_fragmentation_pairwise.xlsx'));

% 2x2 line-plot figure: rows = metrics, cols = detectors (Tampere left, Doorn right)
figF = figure('Color','w','Position',[80 80 1250 860]);
tlF  = tiledlayout(figF, numel(frag_metrics), numel(frag_dets), ...
    'TileSpacing','compact','Padding','compact');
lgdHandles = [];
for im = 1:numel(frag_metrics)
    mname = frag_metrics{im};
    for jd = 1:numel(frag_dets)
        det  = frag_dets{jd};
        Tdet = T(strcmp(T.Detector,det), :);
        ax = nexttile(tlF); hold(ax,'on');
        handles = gobjects(numel(cell_lines),1);
        for ic = 1:numel(cell_lines)
            cl = cell_lines{ic};
            [means_, ci_low, ci_up] = mean_ci_per_div(Tdet, cl, mname, DIVs, alpha, outlierK);
            valid = ~isnan(means_);
            if any(valid)
                xfill = [DIVs(valid), fliplr(DIVs(valid))];
                yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
                fill(ax, xfill, yfill, palette(ic,:), 'FaceAlpha',0.15, ...
                    'EdgeColor','none','HandleVisibility','off');
            end
            handles(ic) = plot(ax, DIVs, means_, '-o', 'Color',palette(ic,:), ...
                'LineWidth',1.8, 'Marker','.','MarkerSize',14,'MarkerFaceColor',palette(ic,:));
        end
        draw_overall_brackets_inax(ax, ovp_det{jd}.(mname), cell_lines, palette, DIVs);
        set(ax,'FontSize',14,'Box','on','XTick',DIVs,'XTickLabel',DIVlabels);
        grid(ax,'on'); ax.GridAlpha = 0.15;
        if im == 1, title(ax, frag_dtitle{jd}, 'FontSize',13,'FontWeight','bold'); end
        if jd == 1, ylabel(ax, frag_mlabel{im}, 'FontSize',15,'FontWeight','bold'); end
        if im == numel(frag_metrics), xlabel(ax,'DIV','FontSize',15); end
        if im == 1 && jd == 1, lgdHandles = handles; end
    end
end
if ~isempty(lgdHandles)
    lgd = legend(lgdHandles, cl_labels, 'Box','off','FontSize',10);
    lgd.Layout.Tile = 'east';
end
% sgtitle(tlF, 'Fragmentation: cell-line trajectories (Tampere vs Doorn)', ...
%     'FontSize',14,'FontWeight','bold');
savebaseF = fullfile(fig_path_plots,'A4_fragmentation_2x2_lineplots');
exportgraphics(figF, [savebaseF '.png'],'Resolution',300);
exportgraphics(figF, [savebaseF '.pdf'],'ContentType','vector');
close(figF);

fprintf('\n===== ALL LME ANALYSES COMPLETED =====\n');
fprintf('Results saved in: %s\n', re_output_path);

%% ====================== LOCAL FUNCTIONS ================================

end

function [p, test_name, eff_size, posthoc] = run_overall_lme(T, metric, fixed_factors, random_factor, cell_lines, DIVs)
posthoc = table(); p = NaN; test_name = ''; eff_size = NaN;
try
    Tlme = T;
    Tlme.(fixed_factors{1}) = categorical(Tlme.(fixed_factors{1}));
    Tlme.(fixed_factors{2}) = categorical(Tlme.(fixed_factors{2}));
    Tlme.(random_factor)    = categorical(Tlme.(random_factor));
    formula = sprintf('%s ~ %s * %s + (1|%s)', ...
        metric, fixed_factors{1}, fixed_factors{2}, random_factor);
    lme = fitlme(Tlme, formula);
    res = residuals(lme);
    try; [~,pn] = swtest(res,0.05);
    catch; try; [~,pn] = lillietest(res); catch; pn = 1; end; end
    anv = anova(lme);
    p   = anv.pValue(end);
    test_name = ternary(pn < 0.01, ...
        'LME (residuals non-normal; interpret with caution)', ...
        'LME (interaction p)');
    try; eff_size = lme.Rsquared.Adjusted; catch; eff_size = NaN; end
    try; posthoc = lme_pairwise_posthoc(lme, fixed_factors{1}, cell_lines);
    catch; posthoc = table(); end
    return;
catch ME
    warning('LME failed: %s. Falling back to Kruskal-Wallis.', ME.message);
end
try
    [p,~,~]   = kruskalwallis(T.(metric), T.(fixed_factors{1}), 'off');
    test_name = 'Kruskal-Wallis (fallback)';
    k         = numel(unique(T.(fixed_factors{1}))); n = numel(T.(metric));
    H         = chi2inv(1-p, k-1);
    eff_size  = (H-k+1)/(n-k);
catch
    p = NaN; test_name = 'FAILED'; eff_size = NaN;
end
end

function [p, test_name, d_eff] = run_overall_lme_detector(T, metric, random_factor, lvl1, lvl2)
p = NaN; test_name = ''; d_eff = NaN;
Tsub = T(ismember(T.Detector,{lvl1,lvl2}) & ~isnan(T.(metric)), :);
if height(Tsub) < 6; return; end
try
    Tlme          = Tsub;
    Tlme.Detector = categorical(Tlme.Detector, {lvl1,lvl2});
    Tlme.(random_factor) = categorical(Tlme.(random_factor));
    lme = fitlme(Tlme, sprintf('%s ~ Detector + (1|%s)', metric, random_factor));
    anv = anova(lme);
    p   = anv.pValue(2);
    test_name = 'LME (Detector effect)';
catch ME
    warning('LME (detector) failed: %s. Falling back to Mann-Whitney.', ME.message);
    x1 = Tsub.(metric)(strcmp(Tsub.Detector,lvl1));
    x2 = Tsub.(metric)(strcmp(Tsub.Detector,lvl2));
    try; p = ranksum(x1,x2); test_name = 'Mann-Whitney (fallback)';
    catch; p = NaN; test_name = 'FAILED'; end
end
x1 = Tsub.(metric)(strcmp(Tsub.Detector,lvl1));
x2 = Tsub.(metric)(strcmp(Tsub.Detector,lvl2));
sA = std(x1,'omitnan'); sB = std(x2,'omitnan');
nA = sum(~isnan(x1));   nB = sum(~isnan(x2));
sp = sqrt(((nA-1)*sA^2+(nB-1)*sB^2)/max(nA+nB-2,1));
if sp > 0; d_eff = (mean(x1,'omitnan')-mean(x2,'omitnan'))/sp; end
end

function out = lme_pairwise_posthoc(lme, factor_name, levels)
out = table(); coefs = lme.Coefficients; coefNames = coefs.Name;
pvals = []; ests = []; pairs_lbl = strings(0,1);
nL = numel(levels);
for i = 1:nL
    for j = i+1:nL
        H     = zeros(1, height(coefs));
        iName = sprintf('%s_%s', factor_name, levels{i});
        jName = sprintf('%s_%s', factor_name, levels{j});
        iIdx  = find(strcmp(coefNames, iName),1);
        jIdx  = find(strcmp(coefNames, jName),1);
        if isempty(iIdx) && isempty(jIdx); continue; end
        if isempty(iIdx);     H(jIdx) = -1;
        elseif isempty(jIdx); H(iIdx) =  1;
        else; H(iIdx) = 1; H(jIdx) = -1; end
        try; p_pair = coefTest(lme, H); catch; continue; end
        est = H * coefs.Estimate;
        pvals(end+1,1) = p_pair; %#ok
        ests(end+1,1)  = est;    %#ok
        pairs_lbl(end+1,1) = sprintf('%s_vs_%s', levels{i}, levels{j});
    end
end
if isempty(pvals); return; end
padj = fdr_bh(pvals);
out  = table(pairs_lbl, ests, pvals, padj, ...
    'VariableNames',{'Pair','Estimate_AminusB','p_raw','p_adj_FDR'});
end

function [p, test_name, d_eff] = run_pairwise_test(xA, xB, paired)
xA = xA(~isnan(xA)); xB = xB(~isnan(xB));
p = NaN; test_name = 'insufficient_n'; d_eff = NaN;
if numel(xA) < 3 || numel(xB) < 3; return; end
nA_ok = is_normal(xA); nB_ok = is_normal(xB);
if nA_ok && nB_ok
    try; equal_var = vartest2(xA,xB) > 0.05; catch; equal_var = false; end
    if paired;        [~,p] = ttest(xA,xB);                     test_name = 'paired t-test';
    elseif equal_var; [~,p] = ttest2(xA,xB);                    test_name = 't-test (equal var)';
    else;             [~,p] = ttest2(xA,xB,'Vartype','unequal'); test_name = 'Welch t-test'; end
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

function s = sig_stars(p)
s = strings(numel(p),1);
for i = 1:numel(p)
    if isnan(p(i));    s(i) = "";
    elseif p(i)<0.001; s(i) = "***";
    elseif p(i)<0.01;  s(i) = "**";
    elseif p(i)<0.05;  s(i) = "*";
    else;              s(i) = "ns";
    end
end
end

function [br, bd, rs, synch] = compute_activity_metrics(tbl, SampFreq, ParamsBurst, recDur, synchThr)
% Per-well activity metrics from the spike table and single-electrode bursts:
%   br    : SEB rate  = mean over burst-detecting electrodes of bursts/min
%   bd    : SEB duration = mean burst duration over all bursts [s]
%   rs    : %% random spikes = 100 * (spikes outside bursts / total spikes)
%   synch : network synchrony = number of active electrodes (FR >= synchThr)
br = NaN; bd = NaN; rs = NaN; synch = NaN;
if ~istable(tbl) || height(tbl)==0 || ...
   ~ismember('Time',tbl.Properties.VariableNames) || ...
   ~ismember('Channel',tbl.Properties.VariableNames)
    return
end
Time    = tbl.Time;
Channel = string(tbl.Channel);
valid   = isfinite(Time) & Time >= 0;
Time    = Time(valid); Channel = Channel(valid);
if isempty(Time), return; end

t0      = min(Time);
recMin  = recDur/60;
elecs   = unique(Channel);
nActive = 0;
brPerElec = [];
allDurs   = [];
sumTotal  = 0;
sumInBurst = 0;

for e = 1:numel(elecs)
    spT = sort(Time(Channel == elecs(e)) - t0);
    if (numel(spT)/recDur) >= synchThr, nActive = nActive + 1; end
    sps = unique(round(spT * SampFreq) + 1);
    sps = sps(sps >= 1);
    sumTotal = sumTotal + numel(sps);
    if numel(sps) < 3, continue; end
    bt = DetectBurstsNeuro(sps, SampFreq, recDur, ParamsBurst);
    if isempty(bt) || size(bt,1)==0, continue; end
    brPerElec(end+1) = size(bt,1) / recMin; %#ok<AGROW>
    allDurs = [allDurs; (bt(:,2)-bt(:,1))/SampFreq]; %#ok<AGROW>
    sumInBurst = sumInBurst + sum(bt(:,3));
end

synch = nActive;
if ~isempty(brPerElec), br = mean(brPerElec, 'omitnan'); end
if ~isempty(allDurs),   bd = mean(allDurs,  'omitnan'); end
if sumTotal > 0,        rs = 100 * (1 - sumInBurst/sumTotal); end
end

function v = compute_lvr(I, R)
% compute_lvr  Revised local variation LvR (Shinomoto et al., 2009) of an
%   interval sequence I, with refractoriness constant R (same units as I).
%   Here applied to the inter-NB-interval sequence. NaN if fewer than 2 intervals.
I = I(:); I = I(isfinite(I) & I > 0);
n = numel(I);
if n < 2, v = NaN; return; end
acc = 0;
for i = 1:n-1
    a = I(i); b = I(i+1);
    acc = acc + (1 - 4*a*b/((a+b)^2)) * (1 + 4*R/(a+b));
end
v = 3/(n-1) * acc;
end

function mfr = compute_mfr_like_mfr_burstR_NBR(spikeTable, recording_dur_s)
% MFR aligned to mfr_burstR_NBR_CI_fixed.m
% Definition: total spikes in the well / recording duration (Hz)

mfr = NaN;

if ~istable(spikeTable) || height(spikeTable) == 0
    return;
end

if ~ismember('Time', spikeTable.Properties.VariableNames) || ...
   ~ismember('Channel', spikeTable.Properties.VariableNames)
    return;
end

mfr = height(spikeTable) / recording_dur_s;

if ~isfinite(mfr)
    mfr = NaN;
end
end

function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end

function [means_, ci_low, ci_up] = mean_ci_per_div(Tsub, cl, metric, DIVs, alpha, outlierK)
means_ = nan(numel(DIVs),1);
ci_low = nan(numel(DIVs),1);
ci_up  = nan(numel(DIVs),1);
for id = 1:numel(DIVs)
    x = Tsub.(metric)(strcmp(Tsub.CellLine,cl) & Tsub.DIV==DIVs(id));
    x = x(~isnan(x));
    x = remove_outliers_iqr(x, outlierK);
    if numel(x) >= 2
        m  = mean(x); s = std(x);
        tc = tinv(0.975, numel(x)-1);
        ci = tc * s / sqrt(numel(x));
        means_(id) = m;
        ci_low(id) = max(m - ci, 0);
        ci_up(id)  = m + ci;
    elseif numel(x)==1
        means_(id) = x; ci_low(id) = x; ci_up(id) = x;
    end
end
end

function h = plot_detector_trace(ax, T, cl, det, metric, DIVs, color, alpha, outlierK)
Tsub = T(strcmp(T.CellLine,cl) & strcmp(T.Detector,det), :);
[means_, ci_low, ci_up] = mean_ci_per_div(Tsub, cl, metric, DIVs, alpha, outlierK);
valid = ~isnan(means_);
if any(valid)
    xfill = [DIVs(valid), fliplr(DIVs(valid))];
    yfill = [ci_low(valid)', fliplr(ci_up(valid)')];
    fill(ax, xfill, yfill, color,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end
h = plot(ax, DIVs, means_,'-o', ...
    'Color',color, 'LineWidth',1.8, ...
    'Marker','.', 'MarkerSize',14, 'MarkerFaceColor',color);
end

function clean = remove_outliers_iqr(data, k)
data = data(isfinite(data));
if numel(data) < 4; clean = data; return; end
q1 = prctile(data, 25); q3 = prctile(data, 75);
iqr_ = q3 - q1;
clean = data(data >= q1 - k*iqr_ & data <= q3 + k*iqr_);
end

function out = compute_overall_cellline_pairs(Tsub, metrics, cell_lines)
out   = struct();
pairs = nchoosek(1:numel(cell_lines),2);
for im = 1:numel(metrics)
    mname = metrics{im};
    npair = size(pairs,1);
    ptbl = nan(npair,1); dtbl = nan(npair,1); mdtbl = nan(npair,1);
    for ip = 1:npair
        A  = cell_lines{pairs(ip,1)}; B = cell_lines{pairs(ip,2)};
        xA = Tsub.(mname)(strcmp(Tsub.CellLine,A));
        xB = Tsub.(mname)(strcmp(Tsub.CellLine,B));
        xA = xA(~isnan(xA)); xB = xB(~isnan(xB));
        if numel(xA)<3 || numel(xB)<3; continue; end
        try; [~,pn1] = swtest(xA,0.05); nA = pn1>0.05;
             [~,pn2] = swtest(xB,0.05); nB = pn2>0.05;
        catch; nA = true; nB = true; end
        if nA && nB; [~,p] = ttest2(xA,xB,'Vartype','unequal');
        else; p = ranksum(xA,xB); end
        ptbl(ip) = p;
        sA = std(xA); sB = std(xB); nAn = numel(xA); nBn = numel(xB);
        sp = sqrt(((nAn-1)*sA^2+(nBn-1)*sB^2)/(nAn+nBn-2));
        if sp > 0; dtbl(ip) = (mean(xA)-mean(xB))/sp; end
        mdtbl(ip) = mean(xA)-mean(xB);
    end
    padj      = fdr_bh(ptbl);
    sig_pairs = struct('idxA',{},'idxB',{},'A',{},'B',{}, ...
                       'p_adj',{},'d',{},'direction',{});
    for ip = 1:npair
        if ~isnan(padj(ip)) && padj(ip) < 0.05
            A = cell_lines{pairs(ip,1)}; B = cell_lines{pairs(ip,2)};
            if mdtbl(ip)>0; dirstr=sprintf('%s>%s',A,B); else; dirstr=sprintf('%s<%s',A,B); end
            sig_pairs(end+1).idxA    = pairs(ip,1); %#ok
            sig_pairs(end).idxB      = pairs(ip,2);
            sig_pairs(end).A         = A;
            sig_pairs(end).B         = B;
            sig_pairs(end).p_adj     = padj(ip);
            sig_pairs(end).d         = dtbl(ip);
            sig_pairs(end).direction = dirstr;
        end
    end
    out.(mname) = sig_pairs;
end
end

function draw_overall_brackets_inax(ax, sig_pairs, cell_lines, palette, DIVs)
% Draw overall cell-line pairwise significance brackets in the right margin of
% the same axes (suitable for subplot/tiled layouts).
if isempty(sig_pairs), return; end
lh = findall(ax,'Type','line','-and','LineStyle','-');
keep = false(size(lh)); for ii = 1:numel(lh), keep(ii) = strcmp(lh(ii).Marker,'.'); end
lh = lh(keep); lh = flipud(lh);
nCL = numel(cell_lines);
if numel(lh) < nCL, return; end
yEnd = nan(nCL,1);
for ic = 1:nCL
    yd = lh(ic).YData; ok = ~isnan(yd);
    if any(ok), yEnd(ic) = yd(find(ok,1,'last')); end
end
xL = xlim(ax); spanX = diff(xL); xR = max(DIVs);
n = numel(sig_pairs); col = 0.07*spanX;
xlim(ax,[xL(1), xR + (n+0.6)*col]);
star = @(pp) ternary_local(pp<0.001,'***',ternary_local(pp<0.01,'**','*'));
for k = 1:n
    sp = sig_pairs(k); yA = yEnd(sp.idxA); yB = yEnd(sp.idxB);
    if isnan(yA) || isnan(yB), continue; end
    cA = palette(sp.idxA,:); cB = palette(sp.idxB,:);
    xv = xR + k*col;
    plot(ax,[xR xv],[yA yA],'-','Color',cA,'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,[xR xv],[yB yB],'-','Color',cB,'LineWidth',1.2,'HandleVisibility','off');
    plot(ax,[xv xv],sort([yA yB]),'-','Color',[0.2 0.2 0.2],'LineWidth',1.2,'HandleVisibility','off');
    text(ax,xv,(yA+yB)/2,star(sp.p_adj),'FontSize',11,'FontWeight','bold', ...
        'Color',[0.1 0.1 0.1],'HorizontalAlignment','center','VerticalAlignment','bottom','Clipping','off');
end
end

function draw_overall_brackets_cellline(fig, ax, sig_pairs, cell_lines, palette, DIVs) %#ok
figure(fig);
lineHandles = findall(ax,'Type','line','-and','LineStyle','-');
keep = false(size(lineHandles));
for ii = 1:numel(lineHandles); keep(ii) = strcmp(lineHandles(ii).Marker,'.'); end
lineHandles = lineHandles(keep); lineHandles = flipud(lineHandles);
nCL = numel(cell_lines);
if numel(lineHandles) < nCL
    warning('draw_overall_brackets_cellline: not enough line handles found.'); return;
end
yEnd = nan(nCL,1);
for ic = 1:nCL
    yd = lineHandles(ic).YData; ok = ~isnan(yd);
    if any(ok); idx = find(ok,1,'last'); yEnd(ic) = yd(idx); end
end
ax_pos   = ax.Position;
panel_x0 = ax_pos(1) + ax_pos(3) + 0.005;
axp = axes('Parent',fig,'Position',[panel_x0, ax_pos(2), 0.20, ax_pos(4)]);
hold(axp,'on'); axp.Visible = 'off';
yLim = ylim(ax); ylim(axp,yLim); xlim(axp,[0 1]);
if isempty(sig_pairs)
    text(axp,0.5,mean(yLim),'no sig. pairs', ...
        'HorizontalAlignment','center','FontSize',8,'Color',[0.55 0.55 0.55]);
    return;
end
n        = numel(sig_pairs); tier_w = 1/n;
star_str = @(p) ternary_local(p<0.001,'***',ternary_local(p<0.01,'**','*'));
for k = 1:n
    sp = sig_pairs(k);
    yA = yEnd(sp.idxA); yB = yEnd(sp.idxB);
    if isnan(yA) || isnan(yB); continue; end
    cA = palette(sp.idxA,:); cB = palette(sp.idxB,:);
    x_ls = (k-1)*tier_w + 0.05*tier_w;
    x_v  = (k-1)*tier_w + 0.55*tier_w;
    x_s  = (k-1)*tier_w + 0.65*tier_w;
    plot(axp,[x_ls,x_v],[yA yA],'-','Color',cA,'LineWidth',1.4);
    plot(axp,[x_ls,x_v],[yB yB],'-','Color',cB,'LineWidth',1.4);
    plot(axp,[x_v x_v],sort([yA yB]),'-','Color',[0.2 0.2 0.2],'LineWidth',1.4);
    text(axp,x_s,(yA+yB)/2,star_str(sp.p_adj), ...
        'FontSize',12,'FontWeight','bold','Color',[0.1 0.1 0.1], ...
        'HorizontalAlignment','left','VerticalAlignment','middle');
end
end

function out = ternary_local(cond, a, b)
if cond; out = a; else; out = b; end
end

function draw_overall_bracket_detector(fig, ax, p_adj, lblA, lblB, colA, colB)
figure(fig);
if p_adj<0.001; stars='***'; elseif p_adj<0.01; stars='**'; else; stars='*'; end
lineHandles = findall(ax,'Type','line','-and','LineStyle','-');
keep = false(size(lineHandles));
for ii = 1:numel(lineHandles); keep(ii) = strcmp(lineHandles(ii).Marker,'.'); end
lineHandles = lineHandles(keep); lineHandles = flipud(lineHandles);
if numel(lineHandles) < 2
    warning('draw_overall_bracket_detector: not enough line handles found.'); return;
end
xEnd = nan(2,1); yEnd = nan(2,1);
for ii = 1:2
    xd = lineHandles(ii).XData; yd = lineHandles(ii).YData; ok = ~isnan(yd);
    if any(ok); idx = find(ok,1,'last'); xEnd(ii)=xd(idx); yEnd(ii)=yd(idx); end
end
if any(isnan(yEnd)); return; end
xLim_orig = xlim(ax); yLim_orig = ylim(ax);
spanX = xLim_orig(2)-xLim_orig(1);
bracketArea_x0 = max(xEnd) + 0.04*spanX;
xlim(ax,[xLim_orig(1), bracketArea_x0 + 0.20*spanX]);
ylim(ax,yLim_orig);
x_v = bracketArea_x0 + 0.55*0.20*spanX;
plot(ax,[xEnd(1),x_v],[yEnd(1) yEnd(1)],'-','Color',colA,'LineWidth',1.2,'HandleVisibility','off');
plot(ax,[xEnd(2),x_v],[yEnd(2) yEnd(2)],'-','Color',colB,'LineWidth',1.2,'HandleVisibility','off');
plot(ax,[x_v x_v],sort([yEnd(1) yEnd(2)]),'-','Color',[0.2 0.2 0.2],'LineWidth',1.4,'HandleVisibility','off');
yMid = (yEnd(1)+yEnd(2))/2;
text(ax,x_v+0.04*spanX,yMid,stars, ...
    'FontSize',12,'FontWeight','bold','Color',[0.1 0.1 0.1], ...
    'HorizontalAlignment','left','VerticalAlignment','middle');
text(ax,x_v+0.04*spanX,yMid-0.05*(yLim_orig(2)-yLim_orig(1)), ...
    sprintf('%s vs %s',lblA,lblB), ...
    'FontSize',7,'Color',[0.35 0.35 0.35], ...
    'HorizontalAlignment','left','VerticalAlignment','top');
end

% HEATMAP FUNCTIONS
function plot_pairwise_heatmap(P, metric, cell_lines, DIVs, DIVlabels, ttl, savebase)
pairs      = nchoosek(1:numel(cell_lines),2); npair = size(pairs,1); ndiv = numel(DIVs);
M          = nan(npair,ndiv);
Lstars     = strings(npair,ndiv);
Ldir       = strings(npair,ndiv);
pairLabels = strings(npair,1);
for ip = 1:npair
    A = cell_lines{pairs(ip,1)}; B = cell_lines{pairs(ip,2)};
    pairLabels(ip) = sprintf('%s vs %s',A,B);
    for id = 1:ndiv
        mask = strcmp(P.Metric,metric) & strcmp(P.GroupA,A) & strcmp(P.GroupB,B) & P.DIV==DIVs(id);
        if ~any(mask); continue; end
        pad = P.p_adj(mask); md = P.MeanDiff_AminusB(mask);
        if isnan(pad)||pad>=0.05; M(ip,id)=0;
        else
            M(ip,id)      = min(-log10(pad),6);
            Lstars(ip,id) = get_star_str(pad);
            Ldir(ip,id)   = ternary_local(md>0,sprintf('%s>%s',A,B),sprintf('%s<%s',A,B));
        end
    end
end
plot_logp_heatmap(M, Lstars, Ldir, cellstr(pairLabels), DIVlabels, ...
    ttl, 'DIV', 'Cell line pair', savebase, true);
end

function plot_detector_heatmap(P, metric, cell_lines, DIVs, DIVlabels, ttl, savebase, lblA, lblB)
ncl    = numel(cell_lines); ndiv = numel(DIVs);
M      = nan(ncl,ndiv); Lstars = strings(ncl,ndiv); Ldir = strings(ncl,ndiv);
md_col = P.Properties.VariableNames{find(contains(P.Properties.VariableNames,'MeanDiff_'),1)};
for ic = 1:ncl
    for id = 1:ndiv
        mask = strcmp(P.Metric,metric) & strcmp(P.CellLine,cell_lines{ic}) & P.DIV==DIVs(id);
        if ~any(mask); continue; end
        pad = P.p_adj(mask); md = P.(md_col)(mask);
        if isnan(pad)||pad>=0.05; M(ic,id)=0;
        else
            M(ic,id)      = min(-log10(pad),6);
            Lstars(ic,id) = get_star_str(pad);
            Ldir(ic,id)   = ternary_local(md>0,sprintf('%s>%s',lblA,lblB),sprintf('%s<%s',lblA,lblB));
        end
    end
end
plot_logp_heatmap(M, Lstars, Ldir, cell_lines, DIVlabels, ...
    ttl, 'DIV', 'Cell line', savebase, true);
end

function plot_logp_heatmap(M, Lstars, Ldir, xLabels, yLabels, ttl, xlbl, ylbl, savebase, doSave) %#ok<INUSD>
% Heatmap of -log10(FDR-adjusted p). Cell coding (from the callers):
%   NaN  -> not tested                       (white, blank cell)
%   0    -> tested but non-significant        (lightest band, "ns" label)
%   >0   -> significant, value = -log10(p_adj) capped at 6, stars drawn.
% Graphic style aligned to rerun_script3_heatmaps_thesis (local_draw_heatmap):
%   discrete 7-level blue colormap, fixed caxis 0-6, adaptive-contrast text,
%   grey grid, large side labels, title on, FDR colorbar.
[nrows,ncols] = size(M);
maxVal = 6;
cmap = [ 0.93 0.93 0.93; ...   % 0  (ns / lightest)
         0.86 0.90 0.97; ...   % 1
         0.72 0.81 0.93; ...   % 2
         0.54 0.69 0.86; ...   % 3
         0.34 0.51 0.76; ...   % 4
         0.16 0.32 0.60; ...   % 5
         0.04 0.17 0.41 ];     % >=6 (darkest)

% Boxes a bit smaller, with left room for the (larger) row labels.
figW = max(360, 95*ncols + 320);
figH = max(260, 64*nrows + 170);
fig  = figure('Color','w','Position',[100 100 figW figH]);
ax   = axes('Parent',fig,'Position',[0.27 0.13 0.60 0.75]);
imagesc(ax, M, 'AlphaData', ~isnan(M));   % NaN cells stay transparent (white)
colormap(ax, cmap); clim(ax,[0 maxVal]);
set(ax,'Color','w'); hold(ax,'on');

% Grey grid.
for r = 0.5:1:nrows+0.5; line(ax,[0.5 ncols+0.5],[r r],'Color',[0.75 0.75 0.75],'LineWidth',1.0,'HandleVisibility','off'); end
for c = 0.5:1:ncols+0.5; line(ax,[c c],[0.5 nrows+0.5],'Color',[0.75 0.75 0.75],'LineWidth',1.0,'HandleVisibility','off'); end

% Significance stars in significant cells, "ns" in tested-but-non-sig cells.
% Text colour adapts to the cell colour luminance for contrast.
for r = 1:nrows
    for c = 1:ncols
        v = M(r,c);
        if isnan(v), continue; end
        star = char(Lstars(r,c));
        if ~isempty(star)
            txt = star; fw = 'bold';   fs = 15;
        elseif v == 0
            txt = 'ns';  fw = 'normal'; fs = 12;
        else
            continue;
        end
        frac = min(v,maxVal)/maxVal;
        cIdx = max(1,min(size(cmap,1), round(1+frac*(size(cmap,1)-1))));
        rgb  = cmap(cIdx,:);
        lum  = 0.299*rgb(1) + 0.587*rgb(2) + 0.114*rgb(3);
        if lum < 0.55, tc = [1 1 1]; else, tc = [0 0 0]; end
        text(ax,c,r,txt,'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',fs,'FontWeight',fw,'Color',tc,'Clipping','off','HandleVisibility','off');
    end
end

% Larger side (tick + axis) labels; title on.
set(ax,'YDir','normal','Box','on','LineWidth',1.1,'TickLength',[0 0],'Layer','top', ...
    'XTick',1:ncols,'XTickLabel',yLabels,'YTick',1:nrows,'YTickLabel',xLabels,'FontSize',13);
ax.XAxis.FontSize = 14; ax.XAxis.FontWeight = 'bold';
ax.YAxis.FontSize = 16; ax.YAxis.FontWeight = 'bold';
xlabel(ax,xlbl,'FontSize',14,'FontWeight','bold');
ylabel(ax,ylbl,'FontSize',16,'FontWeight','bold');
title(ax,ttl,'FontSize',14,'FontWeight','bold','Interpreter','none');

cb = colorbar(ax,'eastoutside');
cb.Label.String   = '-log10(FDR-adjusted p)';
cb.Label.FontWeight = 'bold';
cb.FontSize = 11;
cb.Ticks      = 0:1:maxVal;
cb.TickLabels = {'0','1','2','3','4','5',sprintf('>=%d',maxVal)};

if doSave
    exportgraphics(fig,[savebase '.png'],'Resolution',300);
    exportgraphics(fig,[savebase '.pdf'],'ContentType','vector');
    close(fig);
end
end

function s = get_star_str(p)
if p<0.001; s="***"; elseif p<0.01; s="**"; elseif p<0.05; s="*"; else; s=""; end
end

function place_legend_outside(lgd)
% Move a single-axes legend into the top-right figure margin (next to the
% significance brackets) with a small font, so it never overlaps the curves.
% MATLAB's auto-sized width/height are kept, so entries are not stretched.
set(lgd,'Box','off','FontSize',9,'Units','normalized');
drawnow;
p = get(lgd,'Position');     % [x y w h], auto-sized to contents
p(1) = 0.985 - p(3);         % flush against the right edge of the figure
p(2) = 0.965 - p(4);         % flush against the top of the figure
set(lgd,'Position',p);
end
