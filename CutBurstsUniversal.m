function [bstarts_cut, bends_cut, varargout] = CutBurstsUniversal(~, T_start, T_end, Min_Duration, varargin)
% CutBurstsUniversal  Remove bursts shorter than a minimum duration.
%
% INPUT
%   ~            : unused placeholder (kept for app/object-call compatibility)
%   T_start      : burst start times (vector)
%   T_end        : burst end times (vector)
%   Min_Duration : minimum allowed burst duration (same units as T_start/T_end)
%   Channels     : (optional) per-burst channel labels (cell/vector), pruned in sync
%
% OUTPUT
%   bstarts_cut  : start times of the surviving bursts
%   bends_cut    : end times of the surviving bursts
%   varargout{1} : (optional) pruned Channels, if provided as input

if ~any(nargin==4 | nargin==5)
    error('Wrong input! Check number of arguments')
end
if nargin==5
    Channels = varargin{1};
end

bstarts_cut = T_start;
bends_cut   = T_end;

bdurations = bends_cut - bstarts_cut;
MM = find(bdurations < Min_Duration);
bstarts_cut(MM) = [];
bends_cut(MM)   = [];

if nargin == 5
    Channels(MM) = [];
    varargout{1} = Channels;
end
end
