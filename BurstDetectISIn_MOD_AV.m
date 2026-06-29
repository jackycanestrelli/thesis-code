function [Burst, SpikeBurstNumber, dT] = BurstDetectISIn_MOD_AV(Spikes, N, ISI_N)
% BurstDetectISIn_MOD_AV  ISI_N burst detection (Bakkum et al., 2014), modified
%   to also store the contributing channel labels per burst (for "dominating
%   channel" analysis).
% Ref: Bakkum et al., "Parameters for burst detection", Front. Comput.
%      Neurosci. 7:193 (2014), doi:10.3389/fncom.2013.00193.
%
% INPUT
%   Spikes  : struct with fields
%               .Time    spike times [s]
%               .Channel (optional) spike channel labels
%   N       : number of spikes that must fall within ISI_N to start a burst
%   ISI_N   : ISI_N threshold [s]
%
% OUTPUT
%   Burst : struct with fields
%             .T_start burst start times [s]
%             .T_end   burst end times [s]
%             .S       burst size (number of spikes)
%             .C       burst size (number of channels)
%             .CN      contributing channel labels per burst (cell)
%   SpikeBurstNumber : burst index for each spike (-1 if not in a burst)
%   dT               : N-by-nSpikes matrix of forward/backward ISI_N distances

assert(~isempty(ISI_N), 'ISI_N value is empty!')

% Find where the ISI_N burst condition is met (look both directions per spike)
dT = zeros(N,length(Spikes.Time))+inf;
for j = 0:N-1
    dT(j+1,N:length(Spikes.Time)-(N-1)) = Spikes.Time( (N:end-(N-1))+j ) - ...
        Spikes.Time( (1:end-(N-1)*2)+j );
end
Criteria = zeros(size(Spikes.Time));
Criteria( min(dT)<=ISI_N ) = 1;

% Assign burst numbers to each spike
SpikeBurstNumber = zeros(size(Spikes.Time)) - 1;
INBURST = 0;
NUM_    = 0;
NUMBER  = -1;
BL      = 0;
for i = N:length(Spikes.Time)
    if INBURST == 0
        if Criteria(i)
            INBURST = 1;
            NUM_   = NUM_ + 1;
            NUMBER = NUM_;
            BL     = 1;
        end
    else
        if ~Criteria(i)
            INBURST = 0;
            if BL<N
                SpikeBurstNumber(SpikeBurstNumber==NUMBER) = -1;
                NUM_ = NUM_ - 1;
            end
            NUMBER = -1;
        elseif diff(Spikes.Time([i-(N-1) i])) > ISI_N && BL >= N
            % Split consecutive bursts not separated by a tonic spike.
            NUM_   = NUM_ + 1;
            NUMBER = NUM_;
            BL     = 1;
        else
            BL = BL + 1;
        end
    end
    SpikeBurstNumber(i) = NUMBER;
end

% Assemble burst information
MaxBurstNumber = max(SpikeBurstNumber);
Burst.T_start = zeros(1,MaxBurstNumber);
Burst.T_end   = zeros(1,MaxBurstNumber);
Burst.S       = zeros(1,MaxBurstNumber);
Burst.C       = zeros(1,MaxBurstNumber);
Burst.CN      = cell(1,MaxBurstNumber);

for i = 1:MaxBurstNumber
    ID = find( SpikeBurstNumber==i );
    Burst.T_start(i) = Spikes.Time(ID(1));
    Burst.T_end(i)   = Spikes.Time(ID(end));
    Burst.S(i)       = length(ID);
    if isfield( Spikes, 'Channel' )
        Burst.C(i)  = length( unique(Spikes.Channel(ID)) );
        Burst.CN{i} = Spikes.Channel(ID);
    end
end
end
