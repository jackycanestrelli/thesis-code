function [Params, Pdoorn, ParamsBurst] = get_parameter_set(setNum)
% get_parameter_set  Return the detector parameters for one of the four
%   parameter sets used in the thesis. Set 2 is the reference configuration;
%   Sets 1/3/4 are used for the parameter-sensitivity analysis (Chapter 4.4).
%   Values reproduce Tables 3.3 (Tampere NB), 3.4 (Doorn) and 3.5 (single-
%   electrode bursts). Set 1 keeps the historical "buggy" Gaussian kernel
%   (size 6 instead of 6*sigma+1).
%
% INPUT
%   setNum : parameter set index, 1..4 (default 2 = reference)
%
% OUTPUT
%   Params      : Tampere log-ISI-N network-burst detector parameters
%                 (DetectNetworkBurstsNeuro / ...V2)
%   Pdoorn      : Doorn (2024) network-burst detector parameters
%                 (run_NBdetection_doorn / detectNetworkBursts_Doorn2024)
%   ParamsBurst : single-electrode log-ISI burst detector parameters
%                 (DetectBurstsNeuro), in field Log_ISI_Hist
%
% NOTE: sigma and kernel size are expressed in BINS (timeBin = 25 ms), so a
%       Gaussian sigma of 75 ms corresponds to 3 bins.

if nargin < 1 || isempty(setNum), setNum = 2; end
assert(ismember(setNum,1:4), 'setNum must be 1, 2, 3 or 4.');

timeBin = 25e-3;   % s (fixed across all sets)

% ----- Per-set values (columns = Set 1..4) -----
% Tampere network-burst detector (Table 3.3)
T.BinsPerDec  = [20   22   25   30];
T.Span        = [3    3    2    2];
T.sigma_bins  = [3    3    6    10];     % 75/150/250 ms over a 25 ms bin
T.kernelsize  = [6    19   37   61];     % Set 1 = historical buggy kernel
T.MinCount    = [8    10   12   15];
T.MinDuration = [0.35 0.40 0.50 0.60];
T.MergingFrac = [0.20 0.20 0.30 0.40];
T.MinChannels = [3    4    5    6];
T.MinDomChan  = [3    4    5    6];
T.MinChanContr= [0.20 0.25 0.30 0.40];

% Doorn detector (Table 3.4)
D.sigma_bins  = [3    3    6    10];
D.kernelsize  = [6    19   37   61];
D.movWin_s    = [30   30   25   20];

% Single-electrode burst detector (Table 3.5)
B.BinsperDec  = [20   20   25   25];
B.VoidThr     = [0.99 0.99 0.995 0.995];
B.DefaultThr  = [170  170  150  130];    % ms
B.MaxThr      = [1000 1000 800  600];    % ms
B.MinSpikes   = [8    10   12   15];
B.Span        = [5    4    3    3];
B.MPD         = [2    2    2    3];

s = setNum;

% ----- Tampere network-burst detector parameters -----
Params = struct();
Params.LogISINHist.Smoother   = 'lowess';
Params.LogISINHist.BinsPerDec = T.BinsPerDec(s);
Params.LogISINHist.Span       = T.Span(s);
Params.LogISINHist.N          = 20;                 % ISI-N order, fixed across sets
Params.MinCount               = T.MinCount(s);
Params.MinDuration            = T.MinDuration(s);   % s
Params.MergingType            = 'Adaptive';         % fixed across sets
Params.MergingFracOrThres     = T.MergingFrac(s);
Params.MinChannels            = T.MinChannels(s);
Params.MinDomChannels         = T.MinDomChan(s);
Params.MinChannelContribution = T.MinChanContr(s);
Params.timeBin                = timeBin;
Params.kernelsize             = T.kernelsize(s);
Params.sigma                  = T.sigma_bins(s);

% ----- Doorn network-burst detector parameters -----
Pdoorn = struct();
Pdoorn.maxTime          = 600;            % s
Pdoorn.timeBin          = timeBin;        % s
Pdoorn.movingWindow_sec = D.movWin_s(s);  % s
Pdoorn.sigma            = D.sigma_bins(s);
Pdoorn.kernelsize       = D.kernelsize(s);

% ----- Single-electrode burst detector parameters -----
ParamsBurst = struct();
ParamsBurst.Log_ISI_Hist.Smoother   = 'lowess';
ParamsBurst.Log_ISI_Hist.BinsperDec = B.BinsperDec(s);
ParamsBurst.Log_ISI_Hist.VoidThr    = B.VoidThr(s);
ParamsBurst.Log_ISI_Hist.DefaultThr = B.DefaultThr(s);
ParamsBurst.Log_ISI_Hist.MaxThr     = B.MaxThr(s);
ParamsBurst.Log_ISI_Hist.MinSpikes  = B.MinSpikes(s);
ParamsBurst.Log_ISI_Hist.Span       = B.Span(s);
ParamsBurst.Log_ISI_Hist.MPD        = B.MPD(s);
end
