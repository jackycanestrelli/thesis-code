function [bstarts_merged, bends_merged, varargout] = MergeBurstsUniversal(~, T_start, T_end, Merging_Parameter, Mode, varargin)
% MergeBurstsUniversal  Merge consecutive bursts separated by a small gap.
%
% INPUT
%   ~                 : unused placeholder (kept for app/object-call compatibility)
%   T_start           : burst start times (vector)
%   T_end             : burst end times (vector)
%   Merging_Parameter : merge threshold; absolute gap ('Static') or fraction of
%                       burst duration ('Adaptive')
%   Mode              : 'Static' or 'Adaptive'
%   Channels          : (optional) per-burst channel labels (cell), merged in sync
%
% OUTPUT
%   bstarts_merged : start times after merging
%   bends_merged   : end times after merging
%   varargout{1}   : (optional) merged Channels, if provided as input

if ~any(nargin==5 | nargin==6)
    error('Wrong input! Check number of arguments')
end
if nargin==6
    Channels = varargin{1};
end
if strcmp(Mode, 'Adaptive')
    bdurations = T_end - T_start;
end

for i = 2:numel(T_start)
    gap(i-1) = T_start(i) - T_end(i-1);
end

bstarts_merged = T_start;
bends_merged   = T_end;

for ii = 1:numel(gap)
    switch Mode
        case 'Static'
            Cond = gap(ii) <= Merging_Parameter;
        case 'Adaptive'
            Cond = gap(ii) <= Merging_Parameter*bdurations(ii) || ...
                   gap(ii) <= Merging_Parameter*bdurations(ii+1);
    end
    if Cond==true
        bends_merged(ii)     = nan;
        bstarts_merged(ii+1) = nan;
        if nargin==6
            Channels{ii+1} = [Channels{ii}; Channels{ii+1}];
            Channels{ii}   = [];
        end
    end
    clear Cond
end

bends_merged   = bends_merged(~isnan(bends_merged));
bstarts_merged = bstarts_merged(~isnan(bstarts_merged));

if nargin==6
    notempty = ~cellfun(@isempty, Channels);
    varargout{1} = Channels(notempty);
end
end
