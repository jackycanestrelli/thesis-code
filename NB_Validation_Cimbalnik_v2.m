function NB_Validation_Cimbalnik_v2()
% NB_Validation_Cimbalnik_v2  HFO-based NB validation (Cimbalnik framework) with 3-level event stratification.
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

%  NB_Validation_Cimbalnik_v2.m
%  Estensione di NB_Validation_Cimbalnik.m per supportare:
%    (a) Fragments dentro i NB (via nb_sub.fragmentPeakTimes)
%    (b) Stratificazione a 3 livelli degli eventi HFO:
%          - IN_NB & IN_Fragment
%          - IN_NB & OUT_Fragment   (background sincrono ma fuori Fragment)
%          - OUT_NB                 (background di rete, fuori NB)
%    (c) Conteggi aggregati per well su tutti e 3 i livelli, utili per
%        statistiche successive
%  RIFERIMENTI:
%    - Cimbalnik et al. 2018, Ann Clin Transl Neurol 5(9):1062-76
%      (framework HFO detection: bande HG/R/FR, Hilbert+zscore, z>3, >3 cicli)
%    - Hales et al. 2012, Front Neural Circuits 6:29
%      (HFO in colture MEA dissociate, trasferibilità del framework a in vitro)
%  USO: imposta DETECTOR, LINES_TO_PROCESS, FRAG_HALFWIDTH_S in cima, poi run.

%% ===== USER SETTINGS =====
DETECTOR        = 'Tampere';   % 'Tampere' o 'Doorn'
LINES_TO_PROCESS = {'DD1C','DD5A','TUBA','C10902'};

% Mezza larghezza della finestra Fragment (in secondi) costruita attorno al peak time.
% Un Fragment di peak time t_peak diventa l'intervallo [t_peak - FRAG_HALFWIDTH_S,
% t_peak + FRAG_HALFWIDTH_S]. 25 ms e' un buon default per NB di colture iPSC;
% rifai con 15 ms e 50 ms come sensitivity analysis.
FRAG_HALFWIDTH_S = 0.025;

% Mappa cartelle dei raw (modifica come nel tuo script principale)
% Raw per-well recordings (rec_*.mat with well_data, sf). Override RAW_ROOT in base if needed.
if evalin('base', "exist('RAW_ROOT','var')")
    RAW_ROOT = evalin('base','RAW_ROOT');
else
    RAW_ROOT = 'S:\met_narkilahti_neuro_sto-3700\MEA_data_internship\Jacopo\raw';
end
cellLineMap = struct();
cellLineMap.DD1C   = fullfile(RAW_ROOT, 'DD1C');
cellLineMap.DD5A   = fullfile(RAW_ROOT, 'DD5A');
cellLineMap.TUBA   = fullfile(RAW_ROOT, 'TUBA');
cellLineMap.C10902 = fullfile(RAW_ROOT, 'C10902');

% Output
out_folder = fullfile(outRoot, 'all_outputs', 'Validation_Cimbalnik');
if ~exist(out_folder,'dir'), mkdir(out_folder); end

% Bande Cimbalnik (invariate)
bands(1).name = 'HighGamma';  bands(1).range = [65,100];
bands(2).name = 'Ripple';     bands(2).range = [100,250];
bands(3).name = 'FastRipple'; bands(3).range = [250,600];
n_bands = numel(bands);

% Parametri detector eventi (invariati)
zscore_thresh = 3.0;
min_cycles    = 3;
filter_order  = 8;

%% ===== MAIN LOOP =====
ALL_EV = table();   % accumulo eventi di tutte le cell line per i plot finali
for ci = 1:numel(LINES_TO_PROCESS)
    clName  = LINES_TO_PROCESS{ci};
    clDir   = cellLineMap.(clName);

    nb_struct_name = sprintf('NB_%s_%s', DETECTOR, clName);
    if ~evalin('base', sprintf('exist(''%s'',''var'')', nb_struct_name))
        warning('"%s" non in workspace - salto.', nb_struct_name);
        continue
    end
    nb_root = evalin('base', nb_struct_name);

    wellFiles = dir(fullfile(clDir,'rec_*.mat'));
    if isempty(wellFiles)
        warning('Nessun file rec_*.mat in %s', clDir); continue
    end
    
    fprintf('\n========================================\n');
    fprintf('  %s [%s] - %d well files\n', clName, DETECTOR, numel(wellFiles));
    fprintf('  Fragment half-width = %.0f ms\n', FRAG_HALFWIDTH_S*1000);
    fprintf('========================================\n');

    % Accumulatori PER QUESTA LINEA
    % Aggiunti: in_Fragment, Fragment_idx (nel NB), Fragment_global_idx (nel well)
    EVENTS = repmat(struct( ...
        'cellLine',{}, 'DIV',{}, 'plate',{}, 'well',{}, 'band',{}, ...
        'start_s',{}, 'end_s',{}, 'duration_s',{}, ...
        'peak_z',{}, 'peak_t_s',{}, 'peak_freq_Hz',{}, ...
        'mean_amp',{}, 'electrode_idx',{}, ...
        'in_NB',{}, 'NB_idx',{}, ...
        'in_Fragment',{}, 'Fragment_idx_in_NB',{}, 'Fragment_global_idx',{}), 0,1);

    % Aggiunti: n_Fragment_total, time_in_Fragment_s, time_in_NB_out_Frag_s
    WELLSUM = struct('cellLine',{}, 'DIV',{}, 'plate',{}, 'well',{}, ...
        'n_NB',{}, 'n_Fragment_total',{}, ...
        'rec_duration_s',{}, ...
        'time_in_NB_s',{}, 'time_in_Fragment_s',{}, ...
        'time_in_NB_out_Fragment_s',{}, 'time_out_NB_s',{}, ...
        'n_elec',{}, 'sf_Hz',{});

    for fi = 1:numel(wellFiles)
        fname = wellFiles(fi).name;
        tok = regexp(fname,'rec_([^_]+)_([^_]+)_([^.]+)\.mat','tokens','once');
        if isempty(tok), continue; end
        divLabel = tok{1}; plateID = tok{2}; wellID = tok{3};

        % --- check NB esistenti per questo well ---
        try
            nb_sub = nb_root.(plateID).(divLabel).(wellID);
        catch
            fprintf('  [skip NB] %s.%s.%s.%s\n', nb_struct_name,plateID,divLabel,wellID);
            continue
        end
        if ~isfield(nb_sub,'starts') || isempty(nb_sub.starts), continue; end
        NB_on  = nb_sub.starts(:);
        if isfield(nb_sub,'ends'),       NB_off = nb_sub.ends(:);
        elseif isfield(nb_sub,'stops'),  NB_off = nb_sub.stops(:);
        else, fprintf('  [skip] no ends/stops field\n'); continue
        end
        n_NB = numel(NB_on);

        % --- Fragment peak times per NB (cell array, uno per NB) ---
        has_fragments = isfield(nb_sub, 'fragmentPeakTimes') && ...
                        ~isempty(nb_sub.fragmentPeakTimes);
        if has_fragments
            FRAG_PEAKS_CELL = nb_sub.fragmentPeakTimes;
            if numel(FRAG_PEAKS_CELL) ~= n_NB
                warning('  fragmentPeakTimes ha %d celle ma ci sono %d NB - troncamento conservativo', ...
                    numel(FRAG_PEAKS_CELL), n_NB);
            end
        else
            FRAG_PEAKS_CELL = cell(n_NB,1);
        end

        % --- carica raw ---
        fpath = fullfile(clDir,fname);
        try
            tmp = load(fpath,'well_data','sf');
        catch ME
            warning('Errore load %s: %s', fname, ME.message); continue
        end
        wellStruct = tmp.well_data;
        Fs = double(tmp.sf);

        allF = fieldnames(wellStruct);
        elecFields = sort(allF(startsWith(allF,'E')));
        nElec = numel(elecFields);
        if nElec==0, continue; end

        nSamples = numel(wellStruct.(elecFields{1}));
        rec_dur  = nSamples/Fs;

        % --- maschera dentro NB (invariata) ---
        in_NB_mask = false(nSamples,1);
        for k = 1:n_NB
            i1 = max(1,        round(NB_on(k)*Fs));
            i2 = min(nSamples, round(NB_off(k)*Fs));
            in_NB_mask(i1:i2) = true;
        end
        time_in_NB  = sum(in_NB_mask)/Fs;
        time_out_NB = rec_dur - time_in_NB;

        % --- maschera dentro Fragment + tabella di lookup (Fragment global idx) ---
        % Per ogni Fragment del well costruiamo l'intervallo [tp-h, tp+h] e
        % manteniamo: (a) la mask logica per il tagging veloce,
        % (b) un array parallel con: NB di appartenenza, idx nel NB, idx globale.
        in_Fragment_mask = false(nSamples,1);
        % Allocazione conservativa per la lookup table (espandibile)
        FRAG_TABLE = struct('t_start_s',{},'t_end_s',{},'NB_idx',{}, ...
                            'frag_idx_in_NB',{},'global_idx',{});
        n_Frag_total = 0;

        for k = 1:min(n_NB, numel(FRAG_PEAKS_CELL))
            peaks_k = FRAG_PEAKS_CELL{k};
            if isempty(peaks_k), continue; end
            peaks_k = peaks_k(:);
            for fk = 1:numel(peaks_k)
                tp = peaks_k(fk);
                t1 = tp - FRAG_HALFWIDTH_S;
                t2 = tp + FRAG_HALFWIDTH_S;
                i1 = max(1,        round(t1*Fs));
                i2 = min(nSamples, round(t2*Fs));
                if i2 < i1, continue; end
                in_Fragment_mask(i1:i2) = true;

                n_Frag_total = n_Frag_total + 1;
                FRAG_TABLE(end+1).t_start_s     = (i1-1)/Fs;   %#ok<SAGROW>
                FRAG_TABLE(end).t_end_s         = (i2-1)/Fs;
                FRAG_TABLE(end).NB_idx          = k;
                FRAG_TABLE(end).frag_idx_in_NB  = fk;
                FRAG_TABLE(end).global_idx      = n_Frag_total;
            end
        end
        time_in_Frag = sum(in_Fragment_mask)/Fs;

        % "Buco" dentro NB ma fuori Fragment
        in_NB_outFrag_mask = in_NB_mask & ~in_Fragment_mask;
        time_in_NB_outFrag = sum(in_NB_outFrag_mask)/Fs;

        fprintf('  %s|%s|%s  Fs=%.0f  nE=%d  nNB=%d  nFrag=%d  T=%.0fs ... ', ...
            divLabel,plateID,wellID,Fs,nElec,n_NB,n_Frag_total,rec_dur);

        % --- pre-design filtri (uno per banda) ---
        d_filters = cell(n_bands,1);
        for b = 1:n_bands
            flo = bands(b).range(1);
            fhi = min(bands(b).range(2), Fs/2*0.95);
            d_filters{b} = designfilt('bandpassiir', ...
                'FilterOrder',         filter_order, ...
                'HalfPowerFrequency1', flo, ...
                'HalfPowerFrequency2', fhi, ...
                'SampleRate',          Fs, ...
                'DesignMethod',        'butter');
        end

        t_axis = (0:nSamples-1)'/Fs;
        events_this_well = 0;

        % ============ LOOP BANDE ============
        for b = 1:n_bands
            z_max_across_elec = -inf(nSamples,1);
            elec_at_max       = zeros(nSamples,1,'uint8');

            for e = 1:nElec
                sg = wellStruct.(elecFields{e});
                sg = double(sg(:));
                bs = filtfilt(d_filters{b}, sg);
                env = abs(hilbert(bs));
                mu = mean(env);
                sd = std(env);
                if sd < eps, sd = 1; end
                z = (env - mu)/sd;
                upd = z > z_max_across_elec;
                z_max_across_elec(upd) = z(upd);
                elec_at_max(upd)       = uint8(e);
            end

            f_low = bands(b).range(1);
            min_dur_s = min_cycles / f_low;
            min_dur_smp = round(min_dur_s * Fs);

            above = z_max_across_elec > zscore_thresh;
            d_above = diff([0; above; 0]);
            run_starts = find(d_above==1);
            run_ends   = find(d_above==-1) - 1;
            run_lens   = run_ends - run_starts + 1;
            keep = run_lens >= min_dur_smp;
            run_starts = run_starts(keep);
            run_ends   = run_ends(keep);

            for r = 1:numel(run_starts)
                i1 = run_starts(r);
                i2 = run_ends(r);
                seg = z_max_across_elec(i1:i2);
                [pk_z, pk_idx_loc] = max(seg);
                pk_idx = i1 + pk_idx_loc - 1;
                pk_t   = t_axis(pk_idx);

                % --- tagging IN_NB ---
                in_NB_flag = in_NB_mask(pk_idx);
                nb_idx_match = NaN;
                if in_NB_flag
                    hit = find(NB_on <= pk_t & pk_t <= NB_off, 1, 'first');
                    if ~isempty(hit), nb_idx_match = hit; end
                end

                % --- tagging IN_Fragment ---
                % Un evento e' IN_Fragment se il suo peak time cade dentro
                % la finestra estesa di un Fragment. Cerchiamo il Fragment
                % corrispondente per estrarre indici (locale al NB e globale).
                in_Frag_flag = in_Fragment_mask(pk_idx);
                frag_idx_in_NB = NaN;
                frag_global_idx = NaN;
                if in_Frag_flag && ~isempty(FRAG_TABLE)
                    % Match esatto sul Fragment che contiene il peak
                    for ff = 1:numel(FRAG_TABLE)
                        if FRAG_TABLE(ff).t_start_s <= pk_t && ...
                           pk_t <= FRAG_TABLE(ff).t_end_s
                            frag_idx_in_NB  = FRAG_TABLE(ff).frag_idx_in_NB;
                            frag_global_idx = FRAG_TABLE(ff).global_idx;
                            break
                        end
                    end
                end

                pk_freq = mean(bands(b).range);
                e_pk = double(elec_at_max(pk_idx));
                if e_pk >= 1 && e_pk <= nElec
                    sg = double(wellStruct.(elecFields{e_pk})(:));
                    bs = filtfilt(d_filters{b}, sg);
                    env = abs(hilbert(bs));
                    mean_amp = mean(env(i1:i2));
                else
                    mean_amp = NaN;
                end

                EVENTS(end+1).cellLine             = clName;   %#ok<SAGROW>
                EVENTS(end).DIV                    = divLabel;
                EVENTS(end).plate                  = plateID;
                EVENTS(end).well                   = wellID;
                EVENTS(end).band                   = bands(b).name;
                EVENTS(end).start_s                = (i1-1)/Fs;
                EVENTS(end).end_s                  = (i2-1)/Fs;
                EVENTS(end).duration_s             = (i2-i1+1)/Fs;
                EVENTS(end).peak_z                 = pk_z;
                EVENTS(end).peak_t_s               = pk_t;
                EVENTS(end).peak_freq_Hz           = pk_freq;
                EVENTS(end).mean_amp               = mean_amp;
                EVENTS(end).electrode_idx          = e_pk;
                EVENTS(end).in_NB                  = logical(in_NB_flag);
                EVENTS(end).NB_idx                 = nb_idx_match;
                EVENTS(end).in_Fragment            = logical(in_Frag_flag);
                EVENTS(end).Fragment_idx_in_NB     = frag_idx_in_NB;
                EVENTS(end).Fragment_global_idx    = frag_global_idx;

                events_this_well = events_this_well + 1;
            end
        end

        WELLSUM(end+1) = struct( ...   %#ok<SAGROW>
            'cellLine', clName, 'DIV', divLabel, 'plate', plateID, 'well', wellID, ...
            'n_NB', n_NB, 'n_Fragment_total', n_Frag_total, ...
            'rec_duration_s', rec_dur, ...
            'time_in_NB_s', time_in_NB, ...
            'time_in_Fragment_s', time_in_Frag, ...
            'time_in_NB_out_Fragment_s', time_in_NB_outFrag, ...
            'time_out_NB_s', time_out_NB, ...
            'n_elec', nElec, 'sf_Hz', Fs);

        fprintf('events=%d\n', events_this_well);
    end

    if ~isempty(EVENTS)
        T_events = struct2table(EVENTS);
    else
        T_events = table();
    end
    if ~isempty(WELLSUM)
        T_wells = struct2table(WELLSUM);
    else
        T_wells = table();
    end

    save_name = fullfile(out_folder, sprintf('Validation_%s_%s_v2.mat', DETECTOR, clName));
    save(save_name, 'T_events', 'T_wells', 'bands', ...
        'zscore_thresh', 'min_cycles', 'FRAG_HALFWIDTH_S', '-v7.3');
    fprintf('\n>> Salvato: %s  (events=%d, wells=%d)\n', ...
        save_name, height(T_events), height(T_wells));

    if ~isempty(T_events)
        ALL_EV = [ALL_EV; T_events]; %#ok<AGROW>  accumulo per i plot finali
    end
end

fprintf('\n=== VALIDAZIONE v2 COMPLETATA ===\n');
fprintf('Stratificazione a 3 livelli pronta per i plot.\n');

%% ===== PLOT DI SINTESI =====
% Tre soli grafici, costruiti dagli eventi HFO rilevati (tabella ALL_EV):
%   (1) Peak z-score per le 3 bande            -> boxplot
%   (2) CDF di peak_z per OGNI cell line       -> una figura per cell line
%   (3) Ampiezza dell'envelope (mean_amp)      -> boxplot per banda
plot_summary_validation(ALL_EV, DETECTOR, zscore_thresh, out_folder);

end

function plot_summary_validation(ALL_EV, DETECTOR, zscore_thresh, out_folder)
% Genera i tre plot di sintesi richiesti dalla tabella eventi ALL_EV.
if isempty(ALL_EV) || ~ismember('peak_z', ALL_EV.Properties.VariableNames)
    warning('Nessun evento da plottare - plot di sintesi saltati.');
    return
end

plot_folder = fullfile(out_folder, 'plots');
if ~exist(plot_folder,'dir'), mkdir(plot_folder); end

band_order  = {'HighGamma','Ripple','FastRipple'};
band_labels = {'High Gamma','Ripple','Fast Ripple'};
band_colors = [0.85 0.33 0.10; 0.10 0.55 0.85; 0.45 0.20 0.65];
cl_list = unique(cellstr(string(ALL_EV.cellLine)), 'stable');

% ---------- (1) Peak z-score per le 3 bande (boxplot) ----------
gz = []; gpos = [];
for b = 1:numel(band_order)
    zb = ALL_EV.peak_z(strcmp(ALL_EV.band, band_order{b}));
    zb = zb(isfinite(zb));
    gz   = [gz; zb(:)];                 %#ok<AGROW>
    gpos = [gpos; b*ones(numel(zb),1)]; %#ok<AGROW>
end
if ~isempty(gz)
    figZ = figure('Color','w','Position',[80 80 720 520]);
    axZ = axes(figZ); hold(axZ,'on');
    boxplot(axZ, gz, gpos, 'Labels', band_labels, 'Symbol','.');
    yline(axZ, zscore_thresh, '--k', 'z threshold', 'LineWidth',1.2, ...
        'LabelHorizontalAlignment','left');
    ylabel(axZ, 'Peak z-score', 'FontSize',15);
    title(axZ, sprintf('Peak z-score per band — %s', DETECTOR), 'FontSize',13);
    grid(axZ,'on'); box(axZ,'on');
    exportgraphics(figZ, fullfile(plot_folder, sprintf('01_peakZ_byBand_%s.png', DETECTOR)), 'Resolution',300);
    close(figZ);
    fprintf(' Salvato: 01_peakZ_byBand_%s.png\n', DETECTOR);
end

% ---------- (2) CDF di peak_z per OGNI cell line ----------
for ci = 1:numel(cl_list)
    cln = cl_list{ci};
    figC = figure('Color','w','Position',[80 80 720 520]);
    axC = axes(figC); hold(axC,'on');
    hh = gobjects(0); leg = {};
    for b = 1:numel(band_order)
        zb = ALL_EV.peak_z(strcmp(cellstr(string(ALL_EV.cellLine)),cln) & ...
                           strcmp(ALL_EV.band, band_order{b}));
        zb = zb(isfinite(zb));
        if numel(zb) < 2, continue; end
        [f,x] = ecdf(zb);
        hh(end+1) = plot(axC, x, f, '-', 'Color', band_colors(b,:), 'LineWidth',2.2); %#ok<AGROW>
        leg{end+1} = band_labels{b}; %#ok<AGROW>
    end
    if isempty(hh), close(figC); continue; end
    xline(axC, zscore_thresh, '--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    xlabel(axC, 'Peak z-score', 'FontSize',15);
    ylabel(axC, 'Cumulative probability', 'FontSize',15);
    title(axC, sprintf('Peak z-score CDF — %s | %s', cln, DETECTOR), 'FontSize',13);
    legend(axC, hh, leg, 'Location','southeast', 'Box','off');
    ylim(axC,[0 1]); grid(axC,'on'); box(axC,'on');
    exportgraphics(figC, fullfile(plot_folder, sprintf('02_CDF_peakZ_%s_%s.png', cln, DETECTOR)), 'Resolution',300);
    close(figC);
    fprintf(' Salvato: 02_CDF_peakZ_%s_%s.png\n', cln, DETECTOR);
end

% ---------- (3) Ampiezza dell'envelope per banda (boxplot) ----------
if ismember('mean_amp', ALL_EV.Properties.VariableNames)
    ga = []; gap = [];
    for b = 1:numel(band_order)
        ab = ALL_EV.mean_amp(strcmp(ALL_EV.band, band_order{b}));
        ab = ab(isfinite(ab));
        ga  = [ga; ab(:)];                 %#ok<AGROW>
        gap = [gap; b*ones(numel(ab),1)];  %#ok<AGROW>
    end
    if ~isempty(ga)
        figA = figure('Color','w','Position',[80 80 720 520]);
        axA = axes(figA); hold(axA,'on');
        boxplot(axA, ga, gap, 'Labels', band_labels, 'Symbol','.');
        ylabel(axA, 'Envelope amplitude (\muV)', 'FontSize',15);
        title(axA, sprintf('Envelope amplitude per band — %s', DETECTOR), 'FontSize',13);
        grid(axA,'on'); box(axA,'on');
        exportgraphics(figA, fullfile(plot_folder, sprintf('03_envAmp_byBand_%s.png', DETECTOR)), 'Resolution',300);
        close(figA);
        fprintf(' Salvato: 03_envAmp_byBand_%s.png\n', DETECTOR);
    end
end

fprintf('\n>> Plot di sintesi salvati in: %s\n', plot_folder);
end
