function T = decode_validation_mat()
% decode_validation_mat  Tidy the per-event Cimbalnik validation tables into a CSV.
%
% MATLAB port of the former decode_validation_mat.py helper. The Python script
% existed only because numpy/scipy/mat73 cannot decode a MATLAB v7.3 `table`
% object; inside MATLAB the tables load natively, so this function simply reads
% the four Validation_<DETECTOR>_<line>_v2.mat files, assigns each spectral
% event to one of three mutually exclusive classes and writes a tidy CSV.
%
%   IN_FRAG    : event peak inside a fragment           (in_Fragment == true)
%   IN_NB_OUT  : inside a network burst, outside a frag (in_NB & ~in_Fragment)
%   OUT_NB     : outside any network burst              (~in_NB)
%
% INPUT (read from the base workspace, set by the master script):
%   RESULTS_ROOT  output root (default <pwd>/results). The validation .mat
%                 files are read from <RESULTS_ROOT>/all_outputs/Validation_Cimbalnik.
%
% OUTPUT:
%   events_3class_Tampere.csv  (cellLine, DIV, band, class, peak_z, mean_amp)
%   written next to the validation .mat files; also returned as a table T.
%
% If summary_3class_Tampere.csv is present in the same folder, the High-Gamma
% per-class event counts are cross-checked against it and a PASS/FAIL report is
% printed (optional sanity check, mirrors the Python validate() step).

% ----- output root (same convention as the other cleaned scripts) -----
if evalin('base', "exist('RESULTS_ROOT','var')")
    outRoot = evalin('base','RESULTS_ROOT');
else
    outRoot = fullfile(pwd,'results');
end

DETECTOR  = 'Tampere';
LINES     = {'DD1C','DD5A','TUBA','C10902'};
in_folder = fullfile(outRoot, 'all_outputs', 'Validation_Cimbalnik');
BANDS     = {'HighGamma','Ripple','FastRipple'};

% ----- load + classify every event -----
parts = {};
for ci = 1:numel(LINES)
    clName = LINES{ci};
    fpath  = fullfile(in_folder, sprintf('Validation_%s_%s_v2.mat', DETECTOR, clName));
    if ~exist(fpath,'file')
        warning('decode_validation_mat:missing', 'Missing %s - skipped.', fpath);
        continue
    end
    S = load(fpath, 'T_events');
    if ~isfield(S,'T_events') || isempty(S.T_events)
        warning('decode_validation_mat:empty', 'No T_events in %s - skipped.', fpath);
        continue
    end
    ev = S.T_events;

    cls = repmat({'OUT_NB'}, height(ev), 1);
    cls(logical(ev.in_NB) & ~logical(ev.in_Fragment)) = {'IN_NB_OUT'};
    cls(logical(ev.in_Fragment))                      = {'IN_FRAG'};

    parts{end+1} = table( ...
        cellstr(string(ev.cellLine)), cellstr(string(ev.DIV)), ...
        cellstr(string(ev.band)),     cls, ...
        double(ev.peak_z),            double(ev.mean_amp), ...
        'VariableNames', {'cellLine','DIV','band','class','peak_z','mean_amp'}); %#ok<AGROW>
    fprintf('Decoded %-7s events=%5d  bands=%s  DIVs=%d\n', clName, height(ev), ...
        strjoin(unique(cellstr(string(ev.band)))', ','), numel(unique(cellstr(string(ev.DIV)))));
end

if isempty(parts)
    error('decode_validation_mat:noData', ...
        'No validation files found in %s. Run NB_Validation_Cimbalnik_v2.m first.', in_folder);
end
T = vertcat(parts{:});

% ----- optional cross-check vs summary_3class_Tampere.csv (High Gamma) -----
summaryPath = fullfile(in_folder, 'summary_3class_Tampere.csv');
if exist(summaryPath,'file')
    fprintf('\n=== VALIDATION vs summary_3class_Tampere.csv (High Gamma) ===\n');
    s   = readtable(summaryPath);
    ok  = true;
    classes = {'IN_FRAG','IN_NB_OUT','OUT_NB'};
    for ci = 1:numel(LINES)
        line = LINES{ci};
        for k = 1:numel(classes)
            c = classes{k};
            dval = sum(strcmp(T.cellLine,line) & strcmp(T.band,'HighGamma') & strcmp(T.class,c));
            rsel = strcmp(s.cellLine,line) & strcmp(s.class,c);
            rval = sum(s.n_events_HG(rsel));
            flag = 'OK'; if dval ~= rval, flag = 'MISMATCH'; ok = false; end
            fprintf('  %-7s %-10s decoded=%5d summary=%5d  %s\n', line, c, dval, rval, flag);
        end
    end
    if ok, fprintf('\nALL MATCH\n'); else, fprintf('\n*** VALIDATION FAILED ***\n'); end
else
    fprintf('(summary_3class_Tampere.csv not found in %s - skipping cross-check)\n', in_folder);
end

% ----- write tidy CSV -----
outPath = fullfile(in_folder, 'events_3class_Tampere.csv');
writetable(T, outPath);
fprintf('\nSaved %s  rows: %d\n', outPath, height(T));
end
