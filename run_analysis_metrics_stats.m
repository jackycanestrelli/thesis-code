%% run_analysis_metrics_stats.m
% MASTER SCRIPT 2 — ACTIVITY METRICS + ANALYSIS + STATISTICS
% -------------------------------------------------------------------------
% Computes the per-well activity-metric table (CSV) and then calls every
% figure / statistics function in turn, to reproduce all thesis outputs in
% one clean pass.
%
% REQUIRES in the base workspace (produced by run_preprocessing_processing.m,
% or loaded here):
%   DD1C, DD5A, TUBA, C10902                          raw spike data
%   NB_Tampere_<CL>, NB_Tampere_<CL>_V2, NB_Doorn_<CL> NB detection results
%   SampFreq, ParamsBurst                             burst-detection settings
%
% OUTPUT:
%   metrics_per_well.csv/.mat and all figures under RESULTS_ROOT.
%
% NOTE: run this as a SCRIPT so the called figure functions can read the data
%       from the base workspace.
% -------------------------------------------------------------------------

%% ===== CONFIG =====
RESULTS_ROOT = fullfile(pwd, 'results');     % <- edit: output/figure root
cellLines    = {'DD1C','DD5A','TUBA','C10902'};

% Make the output root visible to every figure function (they read it from base)
assignin('base', 'RESULTS_ROOT', RESULTS_ROOT);
% Global figure style: larger, bold axis fonts; titles are off in the figure functions.
set(groot,'defaultAxesFontSize',14,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultAxesLabelFontSizeMultiplier',1.3);
if ~exist(RESULTS_ROOT,'dir'), mkdir(RESULTS_ROOT); end

% Optional: load NB results / data if they are not already in the workspace
% (uncomment and adjust if you run this script on its own)
% NB_DIR = fullfile(RESULTS_ROOT,'NB_cell_line');
% for ic = 1:numel(cellLines)
%     cl = cellLines{ic};
%     load(fullfile(NB_DIR, ['NB_Tampere_' cl '.mat']));     % -> NB_T
%     load(fullfile(NB_DIR, ['NB_Tampere_' cl '_V2.mat']));  % -> NB_T2
%     load(fullfile(NB_DIR, ['NB_Doorn_'   cl '.mat']));     % -> NB_D
%     assignin('base', ['NB_Tampere_' cl],      NB_T);
%     assignin('base', ['NB_Tampere_' cl '_V2'], NB_T2);
%     assignin('base', ['NB_Doorn_'   cl],      NB_D);
% end

%% ===== 1) ACTIVITY METRICS PER WELL -> CSV =====
% Assemble the structs expected by compute_all_metrics_csv from base variables.
cellData = struct(); nbTampere = struct(); nbDoorn = struct();
for ic = 1:numel(cellLines)
    cl = cellLines{ic};
    if evalin('base', sprintf("exist('%s','var')", cl))
        cellData.(cl) = evalin('base', cl);
    end
    if evalin('base', sprintf("exist('NB_Tampere_%s','var')", cl))
        nbTampere.(cl) = evalin('base', ['NB_Tampere_' cl]);
    end
    if evalin('base', sprintf("exist('NB_Doorn_%s','var')", cl))
        nbDoorn.(cl) = evalin('base', ['NB_Doorn_' cl]);
    end
end

cfg = struct('celllines',{cellLines}, 'saveCSV',true, ...
             'saveFolder',fullfile(RESULTS_ROOT,'all_outputs','metrics'), 'csvName','metrics_per_well.csv');
T = compute_all_metrics_csv(cellData, nbTampere, nbDoorn, ...
        evalin('base','SampFreq'), evalin('base','ParamsBurst'), cfg);
assignin('base','metrics_table', T);
fprintf('Metrics table computed: %d rows.\n', height(T));

%% ===== 2) ACTIVITY-METRIC FIGURES (section 1) =====
runStep(@mfr_burstR_NBR);
runStep(@mfr_burstR_NBR_CI_fixed);
runStep(@mean_SEB_duration_CI);
runStep(@mean_perc_rnd_spikes_fixed);
runStep(@plot_synchrony_lineplot_);

%% ===== 3) NETWORK-BURST FIGURES (Tampere vs Doorn / V1 vs V2) =====
runStep(@NB_duration_CV_IBI_CI_fixed);
runStep(@NBrate_boxplot_compare_TvsD);
runStep(@NBrate_boxplot_V1vsV2_2fig);
runStep(@boxplot_TvD_fragcomaprison);
runStep(@lineplots_tutti_def_fixed_);
runStep(@venn_fragments_fixed);
runStep(@ISIhist_threshold_V1vsV2_fix4);

%% ===== 4) RASTERS =====
runStep(@FR_raster_TvD_modified);
runStep(@RASTER_NBandSEB);

%% ===== 5) STATISTICS & SUMMARY FIGURES (read the metrics CSV) =====
runStep(@stats_RNDSPIKES_new);
runStep(@radar_S_inner_pentagon);
runStep(@outputsection4_LME2);   % parameter-set comparison (LME workflow, aligned to script3)
runStep(@script3_statistics_LME_new_MFR_aligned_integrated);

%% ===== 6) NB VALIDATION (Cimbalnik HFO framework) =====
runStep(@NB_Validation_Cimbalnik_v2);
runStep(@NB_Validation_Plots_v4_representative);
runStep(@decode_validation_mat);   % tidy per-event 3-class CSV (events_3class_Tampere.csv)

fprintf('\nAll analyses completed.\n');

%% ===== helper: run one step, isolating failures =====
function runStep(fh)
name = func2str(fh);
fprintf('\n----- %s -----\n', name);
try
    fh();
catch ME
    warning('%s failed: %s', name, ME.message);
end
end
