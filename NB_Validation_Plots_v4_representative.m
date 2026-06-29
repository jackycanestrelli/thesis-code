function NB_Validation_Plots_v4_representative()
% NB_Validation_Plots_v4_representative  Validation plots (3-level stratification) with representative raster selection.
%
% INPUT : required variables are read from the base workspace (set by the
%         master script). Detected dependencies:
%   DD1C, DD5A, TUBA, C10902
% OUTPUT: figure files saved under RESULTS_ROOT (default <pwd>/results).

% Configurable output root for all saved figures (override from the base
% workspace by setting RESULTS_ROOT before calling this function).
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

% Pull required inputs from the base workspace (populated by the master script).
reqVars = {'DD1C','DD5A','TUBA','C10902'};
for reqI = 1:numel(reqVars)
    reqN = reqVars{reqI};
    if evalin('base', sprintf("exist('%s','var')", reqN))
        eval(sprintf('%s = evalin(''base'', ''%s'');', reqN, reqN));
    end
end
clear reqVars reqI reqN
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);

% NB_Validation_Plots_v4_representative.m
% Plot di validazione con stratificazione a 3 livelli:
% IN_FRAG = IN_NB & IN_Fragment (eventi dentro Fragment)
% IN_NB_OUT = IN_NB & ~IN_Fragment (eventi dentro NB ma fuori Fragment)
% OUT_NB = ~IN_NB (eventi fuori NB)
% NOVITA' v4:
% - Fig 6 NON sceglie piu' il well/NB con piu' frammenti
% - Seleziona invece un raster/NB rappresentativo della cell line,
%   coerente con le metriche aggregate mostrate negli altri plot
% - La selezione del well usa distanza dalla mediana delle metriche HG
% - La selezione dell'NB usa distanza dalla mediana intra-well

%% ===== USER SETTINGS =====
DETECTOR = 'Tampere';
LINES_TO_PLOT = {'DD1C','DD5A','TUBA','C10902'};
in_folder = fullfile(outRoot, 'all_outputs', 'Validation_Cimbalnik');
out_folder_plots = fullfile(outRoot, 'all_outputs', 'Validation_Cimbalnik', 'plots');
if ~exist(out_folder_plots,'dir'), mkdir(out_folder_plots); end

% Raw per-well recordings (rec_*.mat with well_data, sf). Override RAW_ROOT in base if needed.
if evalin('base', "exist('RAW_ROOT','var')")
    RAW_ROOT = evalin('base','RAW_ROOT');
else
    RAW_ROOT = 'S:\met_narkilahti_neuro_sto-3700\MEA_data_internship\Jacopo\raw';
end
cellLineMap = struct();
cellLineMap.DD1C = fullfile(RAW_ROOT, 'DD1C');
cellLineMap.DD5A = fullfile(RAW_ROOT, 'DD5A');
cellLineMap.TUBA = fullfile(RAW_ROOT, 'TUBA');
cellLineMap.C10902 = fullfile(RAW_ROOT, 'C10902');

ILLUSTR_PRE_S = 0.2;
ILLUSTR_POST_S = 0.3;
MUA_SIGMA_MS = 25;
MUA_SPIKE_THR_K = 4;
zscore_thresh = 3;

band_names = {'HighGamma','Ripple','FastRipple'};
class_colors = [0.20 0.65 0.30; 0.95 0.65 0.10; 0.55 0.55 0.55];
class_labels = {'IN_FRAG','IN_NB_OUT','OUT_NB'};
line_colors = [0.70 0.10 0.10; 0.85 0.45 0.10; 0.15 0.30 0.70; 0.35 0.60 0.90];

%% ===== LOAD DATA =====
ALL_EV = table();
ALL_WS = table();
for ci = 1:numel(LINES_TO_PLOT)
    clName = LINES_TO_PLOT{ci};
    fpath = fullfile(in_folder, sprintf('Validation_%s_%s_v2.mat', DETECTOR, clName));
    if ~exist(fpath,'file')
        warning('Manca %s - salto.', fpath); continue
    end
    S = load(fpath);
    if isfield(S,'T_events') && ~isempty(S.T_events), ALL_EV = [ALL_EV; S.T_events]; end %#ok
    if isfield(S,'T_wells') && ~isempty(S.T_wells), ALL_WS = [ALL_WS; S.T_wells]; end %#ok
end
if isempty(ALL_EV)
    error('Nessun dato caricato. Esegui prima NB_Validation_Cimbalnik_v2.m');
end

ALL_EV.event_class = repmat({''}, height(ALL_EV), 1);
m_inFrag = ALL_EV.in_Fragment;
m_inNB_outFrag = ALL_EV.in_NB & ~ALL_EV.in_Fragment;
m_outNB = ~ALL_EV.in_NB;
ALL_EV.event_class(m_inFrag) = {'IN_FRAG'};
ALL_EV.event_class(m_inNB_outFrag) = {'IN_NB_OUT'};
ALL_EV.event_class(m_outNB) = {'OUT_NB'};

ALL_WS.key = cellstr(strcat(string(ALL_WS.cellLine),'|',string(ALL_WS.DIV),'|',string(ALL_WS.plate),'|',string(ALL_WS.well)));
ALL_EV.key = cellstr(strcat(string(ALL_EV.cellLine),'|',string(ALL_EV.DIV),'|',string(ALL_EV.plate),'|',string(ALL_EV.well)));

fprintf('Caricati %d eventi totali da %d wells.\n', height(ALL_EV), height(ALL_WS));

%% ===== HG WELL-LEVEL METRICS FOR REPRESENTATIVE SELECTION =====
HG = ALL_EV(strcmp(ALL_EV.band,'HighGamma'),:);
REP = ALL_WS(:, {'cellLine','DIV','plate','well','key','n_NB','n_Fragment_total', ...
    'time_in_Fragment_s','time_in_NB_out_Fragment_s','time_out_NB_s'});
REP.rate_inFrag = nan(height(REP),1);
REP.rate_inNB_out = nan(height(REP),1);
REP.rate_outNB = nan(height(REP),1);
REP.hitNB_pct = nan(height(REP),1);
REP.hitFrag_pct = nan(height(REP),1);
REP.repScore = nan(height(REP),1);

for w = 1:height(REP)
    k = REP.key{w};
    sub_w = HG(strcmp(HG.key,k),:);
    n_iF = sum(strcmp(sub_w.event_class,'IN_FRAG'));
    n_iO = sum(strcmp(sub_w.event_class,'IN_NB_OUT'));
    n_oN = sum(strcmp(sub_w.event_class,'OUT_NB'));

    if ALL_WS.time_in_Fragment_s(w) > 0
        REP.rate_inFrag(w) = n_iF / ALL_WS.time_in_Fragment_s(w);
    end
    if ALL_WS.time_in_NB_out_Fragment_s(w) > 0
        REP.rate_inNB_out(w) = n_iO / ALL_WS.time_in_NB_out_Fragment_s(w);
    end
    if ALL_WS.time_out_NB_s(w) > 0
        REP.rate_outNB(w) = n_oN / ALL_WS.time_out_NB_s(w);
    end

    if ALL_WS.n_NB(w) > 0
        sub_nb = sub_w(sub_w.in_NB & ~isnan(sub_w.NB_idx), :);
        REP.hitNB_pct(w) = 100 * numel(unique(sub_nb.NB_idx)) / ALL_WS.n_NB(w);
    end
    if ALL_WS.n_Fragment_total(w) > 0
        sub_fg = sub_w(sub_w.in_Fragment & ~isnan(sub_w.Fragment_global_idx), :);
        REP.hitFrag_pct(w) = 100 * numel(unique(sub_fg.Fragment_global_idx)) / ALL_WS.n_Fragment_total(w);
    end
end

rep_metric_names = {'n_NB','n_Fragment_total','rate_inFrag','rate_inNB_out','rate_outNB','hitNB_pct','hitFrag_pct'};

for ci = 1:numel(LINES_TO_PLOT)
    clName = LINES_TO_PLOT{ci};
    sel = strcmp(REP.cellLine, clName);
    sub = REP(sel,:);
    if isempty(sub), continue; end

    M = nan(height(sub), numel(rep_metric_names));
    for j = 1:numel(rep_metric_names)
        M(:,j) = sub.(rep_metric_names{j});
    end

    medv = median(M,1,'omitnan');
    madv = mad(M,1,1);
    madv(~isfinite(madv) | madv==0) = 1;
    Z = abs((M - medv) ./ madv);
    score = mean(Z,2,'omitnan');
    REP.repScore(sel) = score;
end

writetable(REP, fullfile(out_folder_plots, sprintf('representative_well_scores_%s.csv', DETECTOR)));

%% ===== FIG 6 — REPRESENTATIVE NB PER CELL LINE =====
fprintf('\nGenerazione Fig 6 rappresentativa...\n');
for ci = 1:numel(LINES_TO_PLOT)
    clName = LINES_TO_PLOT{ci};
    clDir = cellLineMap.(clName);

    ws_cl = REP(strcmp(REP.cellLine, clName),:);
    if isempty(ws_cl)
        fprintf(' [skip Fig6] nessun well per %s\n', clName); continue
    end

    valid_score = isfinite(ws_cl.repScore);
    if any(valid_score)
        ws_pool = ws_cl(valid_score,:);
    else
        ws_pool = ws_cl;
    end

    [~, ord] = sort(ws_pool.repScore, 'ascend', 'MissingPlacement', 'last');
    ws_pool = ws_pool(ord,:);
    best_row = ws_pool(1,:);

    best_div = char(string(best_row.DIV));
    best_plate = char(string(best_row.plate));
    best_well = char(string(best_row.well));

    fprintf(' %s: representative well = %s|%s|%s | score=%.3f\n', ...
        clName, best_div, best_plate, best_well, best_row.repScore);

    nb_struct_name = sprintf('NB_%s_%s', DETECTOR, clName);
    if ~evalin('base', sprintf('exist(''%s'',''var'')', nb_struct_name))
        fprintf(' [skip Fig6] %s non in workspace\n', nb_struct_name); continue
    end

    try
        nb_sub = evalin('base', sprintf('%s.%s.%s.%s', nb_struct_name, best_plate, best_div, best_well));
    catch
        fprintf(' [skip Fig6] percorso NB non trovato\n'); continue
    end
    if ~isfield(nb_sub,'starts') || isempty(nb_sub.starts)
        fprintf(' [skip Fig6] nessun NB nel well selezionato\n'); continue
    end

    NB_on = nb_sub.starts(:);
    if isfield(nb_sub,'ends')
        NB_off = nb_sub.ends(:);
    elseif isfield(nb_sub,'stops')
        NB_off = nb_sub.stops(:);
    else
        fprintf(' [skip Fig6] no ends/stops field\n'); continue
    end
    n_NB_w = numel(NB_on);

    has_frag_field = isfield(nb_sub,'fragmentPeakTimes') && ~isempty(nb_sub.fragmentPeakTimes) && iscell(nb_sub.fragmentPeakTimes);
    n_frags_per_NB = zeros(n_NB_w,1);
    if has_frag_field
        for kk = 1:min(n_NB_w, numel(nb_sub.fragmentPeakTimes))
            if ~isempty(nb_sub.fragmentPeakTimes{kk})
                n_frags_per_NB(kk) = numel(nb_sub.fragmentPeakTimes{kk});
            end
        end
    end
    nb_durs = NB_off - NB_on;

    key_well = sprintf('%s|%s|%s|%s', clName, best_div, best_plate, best_well);
    ev_well = HG(strcmp(HG.key, key_well),:);
    nb_hit_counts = zeros(n_NB_w,1);
    for kk = 1:n_NB_w
        nb_hit_counts(kk) = sum(ev_well.in_NB & ev_well.NB_idx == kk);
    end

    Xnb = [n_frags_per_NB, nb_durs, nb_hit_counts];
    med_nb = median(Xnb,1,'omitnan');
    mad_nb = mad(Xnb,1,1);
    mad_nb(~isfinite(mad_nb) | mad_nb==0) = 1;
    Znb = abs((Xnb - med_nb) ./ mad_nb);
    nb_score = mean(Znb,2,'omitnan');
    [~, best_NB_idx] = min(nb_score);

    frag_peaks_best = [];
    if has_frag_field && best_NB_idx <= numel(nb_sub.fragmentPeakTimes)
        frag_peaks_best = nb_sub.fragmentPeakTimes{best_NB_idx};
        if ~isempty(frag_peaks_best), frag_peaks_best = frag_peaks_best(:); end
    end

    t_nb_on = NB_on(best_NB_idx);
    t_nb_off = NB_off(best_NB_idx);
    fprintf('    representative NB#%d [%.3f-%.3f s] n_frags=%d dur=%.3f s HFOhits=%d\n', ...
        best_NB_idx, t_nb_on, t_nb_off, numel(frag_peaks_best), nb_durs(best_NB_idx), nb_hit_counts(best_NB_idx));

    fname = sprintf('rec_%s_%s_%s.mat', best_div, best_plate, best_well);
    fpath = fullfile(clDir, fname);
    if ~exist(fpath,'file')
        fprintf(' [skip Fig6] file raw non trovato: %s\n', fname); continue
    end

    tmp = load(fpath,'well_data','sf');
    wellStruct = tmp.well_data;
    Fs = double(tmp.sf);
    allF = fieldnames(wellStruct);
    elecFields = sort(allF(startsWith(allF,'E')));
    nElec = numel(elecFields);
    nSamples = numel(wellStruct.(elecFields{1}));

    i_show_s = max(1, round((t_nb_on - ILLUSTR_PRE_S)*Fs));
    i_show_e = min(nSamples, round((t_nb_off + ILLUSTR_POST_S)*Fs));
    n_show = i_show_e - i_show_s + 1;
    t_show = (0:n_show-1)'/Fs;
    t_nb_on_rel = t_nb_on - (i_show_s-1)/Fs;
    t_nb_off_rel = t_nb_off - (i_show_s-1)/Fs;
    frag_rel = frag_peaks_best - (i_show_s-1)/Fs;

    sigma_smp = max(1, round(MUA_SIGMA_MS*1e-3*Fs));
    k_half = 3*sigma_smp;
    k_t = (-k_half:k_half)';
    gauss_k = exp(-k_t.^2/(2*sigma_smp^2));
    gauss_k = gauss_k/sum(gauss_k);
    MUA = zeros(n_show,1);
    active_count = 0;
    for e = 1:nElec
        sg = double(wellStruct.(elecFields{e})(:));
        sn = sg(i_show_s:i_show_e);
        md = median(abs(sg)/0.6745);
        th = MUA_SPIKE_THR_K*md;
        sp = find(sn < -th);
        if isempty(sp), continue; end
        rs = round(0.0015*Fs);
        kp = true(size(sp));
        for si = 2:numel(sp)
            if kp(si-1) && (sp(si)-sp(si-1)) < rs, kp(si) = false; end
        end
        sp = sp(kp);
        st = zeros(n_show,1);
        st(sp) = 1;
        MUA = MUA + conv(st, gauss_k, 'same')*Fs;
        active_count = active_count + 1;
    end

    n_bands_loc = numel(band_names);
    BSIG_SHOW = zeros(n_show, n_bands_loc);
    ENV_SHOW = zeros(n_show, n_bands_loc);
    ZMAT_SHOW = zeros(n_show, n_bands_loc);
    band_ranges_loc = [65,100; 100,250; 250,600];
    for b = 1:n_bands_loc
        flo = band_ranges_loc(b,1);
        fhi = min(band_ranges_loc(b,2), Fs/2*0.95);
        d_f = designfilt('bandpassiir','FilterOrder',8, ...
            'HalfPowerFrequency1',flo,'HalfPowerFrequency2',fhi, ...
            'SampleRate',Fs,'DesignMethod','butter');
        z_max_full = -inf(nSamples,1);
        for e = 1:nElec
            sg = double(wellStruct.(elecFields{e})(:));
            bs = filtfilt(d_f, sg);
            env = abs(hilbert(bs));
            mu = mean(env);
            sd = std(env);
            if sd < eps, sd = 1; end
            z = (env - mu)/sd;
            upd = z > z_max_full;
            z_max_full(upd) = z(upd);
            if e == 1
                BSIG_SHOW(:,b) = bs(i_show_s:i_show_e);
                ENV_SHOW(:,b) = env(i_show_s:i_show_e);
            end
        end
        ZMAT_SHOW(:,b) = z_max_full(i_show_s:i_show_e);
    end

    t_win_s = (i_show_s-1)/Fs;
    t_win_e = (i_show_e-1)/Fs;
    ev_show = ev_well(ev_well.peak_t_s >= t_win_s & ev_well.peak_t_s <= t_win_e,:);
    ev_t_rel = ev_show.peak_t_s - t_win_s;
    ev_class = ev_show.event_class;

    n_panels = n_bands_loc + 2;
    fig6 = figure('Color','w','Position',[60,30,980,820], ...
        'Name',sprintf('Fig6_representative_%s_%s_%s_%s_NB%d', clName, best_div, best_plate, best_well, best_NB_idx));
    tl6 = tiledlayout(n_panels, 1, 'TileSpacing','tight','Padding','compact');

    ax_mua = nexttile(tl6);
    ymax_mua = max(MUA);
    if ~isfinite(ymax_mua) || ymax_mua <= 0, ymax_mua = 1; end
    fill([t_nb_on_rel t_nb_off_rel t_nb_off_rel t_nb_on_rel], [0 0 ymax_mua*1.1 ymax_mua*1.1], [1 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.4); hold on;
    plot(t_show, MUA, 'Color',[0.15,0.35,0.65], 'LineWidth',0.8);
    for fp = 1:numel(frag_rel)
        if frag_rel(fp) >= 0 && frag_rel(fp) <= t_show(end)
            xline(frag_rel(fp), '--', 'Color',[0.1,0.65,0.2], 'LineWidth',1.2);
        end
    end
    ylabel('MUA (Hz)','FontSize',15);
%     title(sprintf('%s | %s | %s | %s | representative NB #%d/%d | n_{frags}=%d | HFOhits=%d | wellScore=%.3f | nbScore=%.3f', ...
%         clName, best_div, best_plate, best_well, best_NB_idx, n_NB_w, numel(frag_peaks_best), nb_hit_counts(best_NB_idx), best_row.repScore, nb_score(best_NB_idx)), ...
%         'FontSize',8,'FontWeight','bold');
    set(ax_mua,'XTickLabel',[],'Box','off');
    xlim([t_show(1), t_show(end)]);

    band_colors_loc = [0.85 0.33 0.10; 0.10 0.55 0.85; 0.45 0.20 0.65];
    ax_bands = gobjects(n_bands_loc,1);
    for b = 1:n_bands_loc
        ax_bands(b) = nexttile(tl6);
        ymin_b = min(BSIG_SHOW(:,b)); ymax_b = max(BSIG_SHOW(:,b));
        if ~isfinite(ymin_b), ymin_b = -1; end
        if ~isfinite(ymax_b), ymax_b = 1; end
        fill([t_nb_on_rel t_nb_off_rel t_nb_off_rel t_nb_on_rel], [ymin_b ymin_b ymax_b ymax_b], [1 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.4); hold on;
        plot(t_show, BSIG_SHOW(:,b), 'Color', band_colors_loc(b,:)*0.7, 'LineWidth',0.6);
        plot(t_show, ENV_SHOW(:,b), 'Color', band_colors_loc(b,:), 'LineWidth',1.4);
        for fp = 1:numel(frag_rel)
            if frag_rel(fp) >= 0 && frag_rel(fp) <= t_show(end)
                xline(frag_rel(fp),'--','Color',[0.1,0.65,0.2],'LineWidth',1.1);
            end
        end
        ylabel(band_names{b},'FontSize',15,'Interpreter','none');
        set(ax_bands(b),'XTickLabel',[],'Box','off');
        xlim([t_show(1), t_show(end)]);
    end

    ax_hz = nexttile(tl6);
    zmax = max(ZMAT_SHOW(:));
    if ~isfinite(zmax), zmax = zscore_thresh*3; end
    imagesc(t_show, 1:n_bands_loc, ZMAT_SHOW', [0, max(zscore_thresh*3, zmax)]); hold on;
    colormap(ax_hz, hot);
    cb = colorbar; cb.Label.String = 'Z-score (max cross-elec)';
    rectangle('Position',[t_nb_on_rel, 0.5, t_nb_off_rel-t_nb_on_rel, n_bands_loc], 'EdgeColor','r','LineWidth',2.2,'LineStyle','-');
    for fp = 1:numel(frag_rel)
        if frag_rel(fp) >= 0 && frag_rel(fp) <= t_show(end)
            plot([frag_rel(fp) frag_rel(fp)], [0.5, n_bands_loc+0.5], '--', 'Color',[0.1,0.85,0.25], 'LineWidth',1.5);
        end
    end

    if ~isempty(ev_t_rel)
        cl_map = containers.Map(class_labels, {class_colors(1,:), class_colors(2,:), class_colors(3,:)});
        band_y_map = containers.Map(band_names, {1,2,3});
        for ev = 1:numel(ev_t_rel)
            cl = ev_class{ev};
            if isKey(cl_map, cl) && isKey(band_y_map, ev_show.band{ev})
                col = cl_map(cl);
                y_pos = band_y_map(ev_show.band{ev});
                plot(ev_t_rel(ev), y_pos, 'o', 'MarkerFaceColor', col, 'MarkerEdgeColor','w', 'MarkerSize',12, 'LineWidth',1.5);
            end
        end
        h_leg = gobjects(numel(class_labels),1);
        for cls = 1:numel(class_labels)
            h_leg(cls) = plot(nan, nan, 'o', 'MarkerFaceColor', class_colors(cls,:), 'MarkerEdgeColor','w', 'MarkerSize',7);
        end
        legend(h_leg, strrep(class_labels,'_','\_'), 'Location','northeast', 'FontSize',7, 'Box','off');
    end

    if max(ZMAT_SHOW(:)) > zscore_thresh
        contour(t_show, 1:n_bands_loc, ZMAT_SHOW', [zscore_thresh zscore_thresh], 'w--', 'LineWidth',1.0);
    end
    set(ax_hz,'YTick',1:n_bands_loc,'YTickLabel',band_names,'FontSize',7,'Box','off');
    xlabel('Time (s)','FontSize',15);
    ylabel('Band','FontSize',15);
%     title('Z-score (max cross-elec)','FontSize',7);
    xlim([t_show(1), t_show(end)]);

    linkaxes([ax_mua; ax_bands; ax_hz],'x');
%     title(tl6, sprintf('Representative illustrative NB — %s | %s | %s | %s | NB #%d', clName, best_div, best_plate, best_well, best_NB_idx), 'FontSize',11, 'FontWeight','bold');

    out6 = fullfile(out_folder_plots, sprintf('06_representative_%s_%s_NB%d_%s.png', clName, DETECTOR, best_NB_idx, best_div));
    exportgraphics(fig6, out6, 'Resolution',300);
    close(fig6);
    fprintf(' Salvata: %s\n', out6);
end

fprintf('\n=== Fig 6 rappresentativa completata ===\n');

%% ===== HELPERS / NOTES =====
% Questo file rigenera la selezione della Fig 6 con criterio rappresentativo.
% Le altre figure del v3 non vengono ricostruite qui.

end
