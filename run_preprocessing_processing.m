%% run_preprocessing_processing.m
% MASTER SCRIPT 1 — PRE-PROCESSING + PROCESSING
% -------------------------------------------------------------------------
% Loads the raw MEA spike data for every cell line, applies the firing-rate
% filter (pre-processing), then runs all three network-burst detectors
% (Tampere V1, Tampere V2, Doorn) for every plate / DIV / well (processing),
% and saves the resulting NB structures.
%
% INPUT  (edit the CONFIG block below):
%   DATA_ROOT    : folder holding <CellLine>.mat files (DD1C.mat, ...), each
%                  containing one nested struct .(plate).(DIV).(well) = table
%                  with columns Time [s] and Channel.
%   RESULTS_ROOT : folder where NB results (and later figures) are written.
%
% OUTPUT:
%   NB_Tampere_<CL>, NB_Tampere_<CL>_V2, NB_Doorn_<CL>  in the workspace and
%   saved as .mat under RESULTS_ROOT/NB_cell_line.
%
% NOTE: run this as a SCRIPT (not a function) so the variables it creates live
%       in the base workspace and can be reused by run_analysis_metrics_stats.m.
% -------------------------------------------------------------------------

%% ===== CONFIG =====
DATA_ROOT    = fullfile(pwd, 'data');       % <- edit: where <CellLine>.mat live
RESULTS_ROOT = fullfile(pwd, 'results');    % <- edit: where to save outputs
cellLines    = {'DD1C','DD5A','TUBA','C10902'};

applyFRfilter = true;                        % pre-processing on/off
SampFreq      = 12500;                       % [Hz] (Axion Maestro, Section 3.1.4)
PARAM_SET     = 2;                           % reference configuration (thesis Set 2)

% Detector parameters for the chosen set (thesis Tables 3.3-3.5):
%   Params      - Tampere log-ISI-N network-burst detector
%   Pdoorn      - Doorn (2024) network-burst detector
%   ParamsBurst - single-electrode log-ISI burst detector (DetectBurstsNeuro)
[Params, Pdoorn, ParamsBurst] = get_parameter_set(PARAM_SET);

saveDir = fullfile(RESULTS_ROOT, 'NB_cell_line');
if ~exist(saveDir,'dir'), mkdir(saveDir); end

%% ===== MAIN LOOP OVER CELL LINES =====
for ic = 1:numel(cellLines)
    cl = cellLines{ic};
    fprintf('\n==================== %s ====================\n', cl);

    % --- Load raw spike data into a variable named after the cell line ---
    f = fullfile(DATA_ROOT, [cl '.mat']);
    if ~exist(f,'file')
        warning('[%s] data file not found: %s — skipping.', cl, f);
        continue
    end
    S = load(f);
    CellData = S.(cl);

    % --- PRE-PROCESSING: firing-rate filter ---
    if applyFRfilter
        CellData = filter_electrodes_by_FR(CellData);
    end
    assignin('base', cl, CellData);   % keep filtered data in base workspace

    % --- PROCESSING: network-burst detection ---
    NB_T   = run_NBdetection_tampere(CellData, Params, @DetectNetworkBurstsNeuro);
    NB_T2  = run_NBdetection_tampere(CellData, Params, @DetectNetworkBurstsNeuroV2);
    NB_D   = run_NBdetection_doorn(CellData, Pdoorn);

    % --- Store with the canonical names used by the analysis scripts ---
    assignin('base', ['NB_Tampere_' cl],        NB_T);
    assignin('base', ['NB_Tampere_' cl '_V2'],   NB_T2);
    assignin('base', ['NB_Doorn_'   cl],         NB_D);

    % --- Save ---
    save(fullfile(saveDir, ['NB_Tampere_' cl '.mat']),      'NB_T',  '-v7.3');
    save(fullfile(saveDir, ['NB_Tampere_' cl '_V2.mat']),   'NB_T2', '-v7.3');
    save(fullfile(saveDir, ['NB_Doorn_'   cl '.mat']),      'NB_D',  '-v7.3');
    fprintf('[%s] detection done and saved.\n', cl);
end

% Expose parameters needed by the analysis stage
assignin('base', 'SampFreq',    SampFreq);
assignin('base', 'ParamsBurst', ParamsBurst);   % single-electrode burst params
assignin('base', 'RESULTS_ROOT', RESULTS_ROOT);

fprintf('\nAll pre-processing and processing completed.\n');
