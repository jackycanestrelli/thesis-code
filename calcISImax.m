function [ISImax, pks, flags] = calcISImax(bins, ISIhist, mpd, voidParamTh, ISITh)
% calcISImax  Find the ISI value that optimally separates intra-burst from
%             inter-burst ISIs on a logarithmic ISI histogram.
% Original author: Valentina Pasquale, Italian Institute of Technology (2008).
% Ref: Pasquale, Martinoia, Chiappalone, J. Comput. Neurosci. (2010),
%      doi:10.1007/s10827-009-0175-1.
%
% INPUT
%   bins        : x values of the logarithmic ISI histogram (vector)
%   ISIhist     : y values of the ISI histogram (vector, same size as bins)
%   mpd         : minimum distance between two peaks, in points (positive integer)
%   voidParamTh : minimum void parameter for two peaks to count as separated, in (0,1)
%   ISITh       : max time window to search the intra-burst peak [ms] (positive integer)
%
% OUTPUT
%   ISImax : ISI separating intra-burst and out-burst ISIs (scalar); -1 if no
%            intra-burst peak below ISITh was found
%   pks    : n-by-2 matrix with [x y] coordinates of detected histogram peaks
%   flags  : 1-by-2 vector; flags(1)=1 if channel is likely bursting,
%            flags(2)=1 if two peaks separated by a valid minimum were found

ISImax = [];
flags  = zeros(1,2);
pks    = [];

% Argument checks
if ~isvector(bins)
    error('The first argument must be a 1-by-N or N-by-1 vector where N >= 1.')
else
    bins = bins(:);
end
if ~isvector(ISIhist)
    error('The second argument must be a 1-by-N or N-by-1 vector where N >= 1.')
else
    ISIhist = ISIhist(:);
end
if length(bins) ~= length(ISIhist)
    error('The first two arguments must be two vectors of the same size.')
end
if ~(isscalar(mpd) && ~(mod(mpd,1))) || mpd <= 0
    error('The third argument must be a positive scalar integer.')
end
if ~isscalar(voidParamTh) || voidParamTh <= 0 || voidParamTh >= 1
    error('The fourth argument must be a scalar value in the range [0,1].')
end
if ISITh <= 0 || ~(isscalar(ISITh) && ~(mod(ISITh,1)))
    error('The fifth argument must be a positive scalar integer.')
end

% Work on the log-x histogram (x already linearly spaced)
xx = bins;
yy = ISIhist;

if (~isempty(yy))
    [peaks,locs] = findpeaks(yy,'MinPeakDistance',mpd);
    if ~isempty(peaks) && any(peaks)
        pks = [xx(locs) peaks(:)];
        numPeaks = size(pks,1);

        % Candidate intra-burst peaks (x < ISITh); keep the highest one
        idxPeakIntraBurst = find(pks(:,1)<ISITh);
        if(numel(idxPeakIntraBurst)>1)
            [~,idxMax] = max(pks(idxPeakIntraBurst,2));
            idxPeakIntraBurst = idxPeakIntraBurst(idxMax);
        else if(isempty(idxPeakIntraBurst))
                ISImax = -1;
                flags(1) = 1;
                flags(2) = 0;
                return;
            end
        end

        y1 = pks(idxPeakIntraBurst,2);
        x1 = pks(idxPeakIntraBurst,1);
        locs1 = find(xx==x1);
        numPeaksAfterBurst = numPeaks-idxPeakIntraBurst;

        if numPeaksAfterBurst == 0
            flags(1) = 1;
            flags(2) = 0;
            return;
        end

        if numPeaksAfterBurst >= 1
            yMin          = zeros(numPeaksAfterBurst-1,1);
            idxMin        = zeros(numPeaksAfterBurst-1,1);
            voidParameter = zeros(numPeaksAfterBurst-1,1);
            c = 0;
            for j = idxPeakIntraBurst:numPeaks
                c = c+1;
                x2 = pks(j,1);
                locs2 = find(xx==x2);
                y2 = pks(j,2);
                [yMin(c),tempIdxMin] = min(yy(locs1:locs2));
                idxMin(c) = tempIdxMin+locs1-1;
                % Void parameter: degree of separation between the two peaks
                voidParameter(c) = 1-(yMin(c)/sqrt(y1.*y2));
            end
            idxMaxVoidParameter = find(voidParameter>=voidParamTh,1);
            if isempty(idxMaxVoidParameter)
                flags(1) = 1;
                flags(2) = 0;
                return;
            end
            ISImax = xx(idxMin(idxMaxVoidParameter));
            flags = [1 1];
        end
    end
end
end
