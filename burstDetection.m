function [burstTrain, burstEventTrain] = burstDetection(peakTrain, ...
    ISImax, flags, ISImaxTh, maxISImaxTh, sf, minNumSpikes)
% burstDetection  Detect bursts on a single spike train (string method).
% Original author: Valentina Pasquale, Italian Institute of Technology (2008).
% Ref: Pasquale, Martinoia, Chiappalone, J. Comput. Neurosci. (2010),
%      doi:10.1007/s10827-009-0175-1.
%
% INPUT
%   peakTrain    : [nSamples x 1] array, 1 where a spike occurs, 0 elsewhere
%   ISImax       : ISI separating intra/out-burst ISIs [ms] (from calcISImax)
%   flags        : 1-by-2 flags vector (from calcISImax)
%   ISImaxTh     : ISI [ms] that switches between single- and double-threshold
%                  detection (core threshold = ISImaxTh, extension = ISImax)
%   maxISImaxTh  : maximum allowed ISImax [ms]; above it ISImaxTh is used
%   sf           : sampling frequency [samples/s]
%   minNumSpikes : minimum number of spikes per burst
%
% OUTPUT
%   burstTrain      : [nBursts x 6+] matrix:
%                     1 start [sample], 2 end [sample], 3 nSpikes,
%                     4 duration [s], 5 inter-burst interval [s],
%                     6 burst period [s], 7 ISI threshold [sample]
%   burstEventTrain : sparse [nSamples x 1], 1 at the first spike of each burst

burstTrain = [];
burstEventTrain = [];

% Input checks
if ~isvector(peakTrain)
    error('The first argument must be a 1-by-N or N-by-1 vector where N >= 1.')
else
    peakTrain = peakTrain(:);
end
if ~all(size(flags) == [1,2]) || any(flags < 0)
    error('The third argument must be a vector of size (1,2) of positive values.')
end
if ~isscalar(ISImaxTh) || ISImaxTh <= 0
    error('The fourth argument must be a positive scalar.')
end
if ~isscalar(maxISImaxTh) || maxISImaxTh <= 0
    error('The fifth argument must be a positive scalar.')
end
if ~(isscalar(sf) && ~(mod(sf,1))) || sf <= 0
    error('The sixth argument must be a positive scalar integer.')
end
if ~(isscalar(minNumSpikes) && ~(mod(minNumSpikes,1))) || minNumSpikes <= 0
    error('The seventh argument must be a positive scalar integer.')
end

% Select single- or double-threshold mode
if flags(2) && ISImax < maxISImaxTh
    if ISImax > ISImaxTh
        ISImaxsample = round(ISImaxTh/1000*sf);
        addWinsample = round(ISImax/1000*sf);
        amplFlag = 1;
    else
        ISImaxsample = round(ISImax/1000*sf);
        amplFlag = 0;
    end
elseif ISImax < 0
    burstTrain = [];
    burstEventTrain = [];
    return
else
    ISImaxsample = round(ISImaxTh/1000*sf);
    amplFlag = 0;
end

timestamp = find(peakTrain);
if ~isempty(timestamp)
    allisi2 = [0;diff(timestamp) <= ISImaxsample;0];   % ISI <= ISImax [samples]
    edges2  = diff(allisi2);                            % burst edges (begin/end)
    edgeup2   = find(edges2 == 1);
    edgedown2 = find(edges2 == -1);

    if ((length(edgedown2)>=2) && (length(edgeup2)>=2))
        if amplFlag
            tsEdgeup2   = timestamp(edgeup2);
            tsEdgedown2 = timestamp(edgedown2);
            tempIBI     = tsEdgeup2(2:end)-tsEdgedown2(1:end-1);
            burst2join  = tempIBI<=addWinsample;
            if any(burst2join)
                edgeup2(find(burst2join)+1)=[];
                edgedown2(burst2join)=[];
            end
            numSpikesInBurst = (edgedown2-edgeup2+1);
            validBursts = numSpikesInBurst >= minNumSpikes;
            if ~isempty(validBursts)
                edgeup2 = edgeup2(validBursts);
                edgedown2 = edgedown2(validBursts);
            end

            allisi1 = [0;diff(timestamp)<=addWinsample;0];
            edges1  = diff(allisi1);
            edgeup1   = find(edges1 == 1);
            edgedown1 = find(edges1 == -1);
            numSpikesInBurst1 = (edgedown1-edgeup1+1);
            validBursts1 = numSpikesInBurst1 >= minNumSpikes;
            if ~isempty(validBursts1)
                edgeup1 = edgeup1(validBursts1);
                edgedown1 = edgedown1(validBursts1);
            end

            allEdgeUp1   = [timestamp(edgeup1) ones(length(edgeup1),1) 1*ones(length(edgeup1),1)];
            allEdgeDown1 = [timestamp(edgedown1) -1*ones(length(edgedown1),1) 1*ones(length(edgeup1),1)];
            allEdge1     = [allEdgeUp1;allEdgeDown1];
            allEdgeUp2   = [timestamp(edgeup2) ones(length(edgeup2),1) 2*ones(length(edgeup2),1)];
            allEdgeDown2 = [timestamp(edgedown2) -1*ones(length(edgedown2),1) 2*ones(length(edgeup2),1)];
            allEdge2     = [allEdgeUp2;allEdgeDown2];
            allEdge      = [allEdge1; allEdge2];
            allEdgeSort  = sortrows(allEdge,1);
            % Columns: [timestamp, edge type (1 rise/-1 fall), train (1 large/2 strict ISI)]

            burstBegin = find(allEdgeSort(:,2)==1 & allEdgeSort(:,3)==1);
            b = 0;
            newBurstDetection = [];
            for ii = 1:length(burstBegin)
                if (allEdgeSort(burstBegin(ii)+1,2)==-1 && allEdgeSort(burstBegin(ii)+1,3)==1)
                    continue
                else
                    thisBurstEnd = find(allEdgeSort(burstBegin(ii):end,2)==-1 & allEdgeSort(burstBegin(ii):end,3)==1,1);
                    thisBurstEnd = thisBurstEnd+burstBegin(ii)-1;
                    subBurstBegin = find(allEdgeSort(burstBegin(ii):thisBurstEnd,2)==1 & allEdgeSort(burstBegin(ii):thisBurstEnd,3)==2);
                    subBurstBegin = subBurstBegin+burstBegin(ii)-1;
                    if length(subBurstBegin)>1
                        subBurstBegin = [burstBegin(ii);subBurstBegin(2:end)];
                        subBurstEnd   = [subBurstBegin(2:end);thisBurstEnd];
                        for jj = 1:length(subBurstBegin)
                            b = b+1;
                            timestampBegin   = allEdgeSort(subBurstBegin(jj),1);
                            tempTimestampEnd = allEdgeSort(subBurstEnd(jj),1);
                            if jj ~= length(subBurstBegin)
                                spks = find(peakTrain(timestampBegin:tempTimestampEnd-1));
                            else
                                spks = find(peakTrain(timestampBegin:tempTimestampEnd));
                            end
                            timestampEnd = timestampBegin+spks(end)-1;
                            numSpikesThisBurst = length(spks);
                            duration_s = (timestampEnd-timestampBegin)./sf;
                            newBurstDetection(b,:) = [timestampBegin timestampEnd numSpikesThisBurst duration_s];
                        end
                    else
                        b = b+1;
                        timestampBegin = allEdgeSort(burstBegin(ii),1);
                        timestampEnd   = allEdgeSort(thisBurstEnd,1);
                        numSpikesThisBurst = sum(spones(peakTrain(timestampBegin:timestampEnd)));
                        duration_s = (timestampEnd-timestampBegin)./sf;
                        newBurstDetection(b,:) = [timestampBegin timestampEnd numSpikesThisBurst duration_s];
                    end
                end
            end
        else
            numSpikesInBurst = (edgedown2-edgeup2+1);
            validBursts = numSpikesInBurst >= minNumSpikes;
            if ~isempty(validBursts)
                edgeup2 = edgeup2(validBursts);
                edgedown2 = edgedown2(validBursts);
            end
            tsEdgeup2   = timestamp(edgeup2);
            tsEdgedown2 = timestamp(edgedown2);
            newBurstDetection = [tsEdgeup2, tsEdgedown2, (edgedown2-edgeup2+1), ...
                (tsEdgedown2-tsEdgeup2)/sf];
        end

        if ~isempty(newBurstDetection)
            binit = newBurstDetection(:,1);
            burstEventTrain = sparse(binit, ones(length(binit),1), 1);
            bp  = [diff(binit)/sf; 0];                                                  % burst period [s]
            ibi = [((newBurstDetection(2:end,1)- newBurstDetection(1:end-1,2))/sf); 0]; % inter-burst interval [s]
            ISIsample = repmat(ISImaxsample,size(newBurstDetection,1),1);
            burstTrain = [newBurstDetection, ibi, bp, ISIsample];
        end
    end
end
end
