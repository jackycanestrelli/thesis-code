function Spikes = filter_electrodes_by_FR(Spikes, varargin)
% filter_electrodes_by_FR  Three-level firing-rate filter for MEA spike data:
%   1) electrode-level : drop electrodes with FR < frThreshold
%   2) well-level (FR) : drop wells with mean FR < wellFRThreshold
%   3) well-level (N)  : drop wells with < minActiveElectrodes active electrodes
%
% INPUT
%   Spikes              : nested struct Spikes.(plate).(div).(well) = table with
%                         columns Time [s] and Channel
%   frThreshold         : (optional) min FR per electrode [Hz],            default 0.2
%   recDuration         : (optional) recording duration [s],              default auto
%   wellFRThreshold     : (optional) min well mean FR [Hz], 0 disables,   default 0.2
%   minActiveElectrodes : (optional) min active electrodes/well, 0 disables, default 3
%
% OUTPUT
%   Spikes : same struct with inactive electrodes removed and inactive wells emptied

frThreshold         = 0.2;
recDuration         = [];
wellFRThreshold     = 0.2;
minActiveElectrodes = 3;

if nargin >= 2 && ~isempty(varargin{1}), frThreshold         = varargin{1}; end
if nargin >= 3 && ~isempty(varargin{2}), recDuration         = varargin{2}; end
if nargin >= 4 && ~isempty(varargin{3}), wellFRThreshold     = varargin{3}; end
if nargin >= 5 && ~isempty(varargin{4}), minActiveElectrodes = varargin{4}; end

fprintf('\n=== Firing Rate Filter ===\n');
fprintf('  Electrode threshold   : %.3f Hz\n', frThreshold);
fprintf('  Well mean FR threshold: %.3f Hz\n', wellFRThreshold);
fprintf('  Min active electrodes : %d\n',      minActiveElectrodes);

totalElecRemoved = 0;
totalElecKept    = 0;
totalWellRemoved = 0;
totalWellKept    = 0;

plateNames = fieldnames(Spikes);

for p = 1:numel(plateNames)
    plate    = plateNames{p};
    divNames = fieldnames(Spikes.(plate));

    for d = 1:numel(divNames)
        div       = divNames{d};
        wellNames = fieldnames(Spikes.(plate).(div));

        for w = 1:numel(wellNames)
            wl  = wellNames{w};
            tbl = Spikes.(plate).(div).(wl);

            if isempty(tbl) || height(tbl) == 0, continue; end

            if isempty(recDuration)
                dur = max(tbl.Time) - min(tbl.Time);
                if dur <= 0, dur = 1; end
            else
                dur = recDuration;
            end

            % Step 1: electrode-level filter
            electrodes = unique(tbl.Channel);
            keepElec   = strings(0);
            elecFRs    = [];

            for e = 1:numel(electrodes)
                elec    = electrodes(e);
                nSpikes = sum(tbl.Channel == elec);
                fr      = nSpikes / dur;

                if fr >= frThreshold
                    keepElec(end+1) = elec; %#ok
                    elecFRs(end+1)  = fr;   %#ok
                else
                    fprintf('  [ELEC REMOVED] %s|%s|%s elec %s — FR=%.4f Hz (%d spikes/%.0fs)\n', ...
                        plate, div, wl, elec, fr, nSpikes, dur);
                    totalElecRemoved = totalElecRemoved + 1;
                end
            end

            totalElecKept = totalElecKept + numel(keepElec);

            if isempty(keepElec)
                Spikes.(plate).(div).(wl) = tbl(false, :);
            else
                mask = ismember(tbl.Channel, keepElec);
                Spikes.(plate).(div).(wl) = tbl(mask, :);
            end

            % Steps 2 + 3: well-level filters
            nActive    = numel(keepElec);
            wellMeanFR = 0;
            if ~isempty(elecFRs), wellMeanFR = mean(elecFRs); end

            failFR = wellFRThreshold     > 0 && wellMeanFR < wellFRThreshold;
            failN  = minActiveElectrodes > 0 && nActive   < minActiveElectrodes;

            if failFR || failN
                reason = '';
                if failFR
                    reason = sprintf('mean FR=%.4f Hz < %.3f Hz', wellMeanFR, wellFRThreshold);
                end
                if failN
                    if ~isempty(reason), reason = [reason ' | ']; end %#ok
                    reason = [reason sprintf('%d active elec < %d required', nActive, minActiveElectrodes)]; %#ok
                end
                fprintf('  [WELL REMOVED]  %s|%s|%s — %s\n', plate, div, wl, reason);
                Spikes.(plate).(div).(wl) = tbl(false, :);
                totalWellRemoved = totalWellRemoved + 1;
            else
                fprintf('  [WELL KEPT]     %s|%s|%s — mean FR=%.4f Hz | active elec=%d\n', ...
                    plate, div, wl, wellMeanFR, nActive);
                totalWellKept = totalWellKept + 1;
            end

        end
    end
end

fprintf('\n--- Summary ---\n');
fprintf('  Electrodes kept   : %d\n', totalElecKept);
fprintf('  Electrodes removed: %d (FR < %.3f Hz)\n', totalElecRemoved, frThreshold);
fprintf('  Wells kept        : %d\n', totalWellKept);
fprintf('  Wells removed     : %d\n', totalWellRemoved);
fprintf('=================================\n\n');
end
