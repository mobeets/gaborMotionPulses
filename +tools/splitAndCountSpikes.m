function [Y, lbins, rbins] = splitAndCountSpikes(spikeTimes, ...
    splitEvent, alignEvent, tL, tR, tWidth, tShift, inds)
% 
    if nargin < 8
        inds = true(numel(splitEvent),1);
    end
    if nargin < 7
        tShift = nan;
    end
    if nargin < 6
        tWidth = nan;
    end
    inds = inds & ~isnan(splitEvent) & ~isnan(alignEvent);
    alignEvent = alignEvent(inds);
    splitEvent = splitEvent(inds);

    t0 = alignEvent - tL;
    t1 = alignEvent + tR;
    [Y, lbins, rbins] = tools.countSpikesWithinWindowByEvent(...
        spikeTimes, splitEvent, t0, t1, tWidth, tShift);
    lbins = lbins - tL;
    rbins = rbins - tL;
end
