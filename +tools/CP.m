function [cp, lbins, Y] = CP(spikeTimes, splitEvent, alignEvent, ...
    tL, tR, tWidth, tShift, inds)
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
    assert(numel(unique(splitEvent))==2);

    t0 = alignEvent - tL;
    t1 = alignEvent + tR;
    [Y, lbins] = countSpikesWithinWindowByEvent(spikeTimes, splitEvent, ...
        t0, t1, tWidth, tShift);
    
    nbins = size(Y{1},2);
    cp = nan(nbins,1);
    for ii = 1:nbins
        cp(ii) = tools.AUC(Y{1}(:,ii), Y{2}(:,ii));
    end
end

function [Y, lbins] = countSpikesWithinWindowByEvent(sps, ev, t0s, t1s, ...
    tWidth, tShift)
% 
    categs = sort(unique(ev));
    Y = cell(numel(categs), 1);
    for ii = 1:numel(categs)
        inds = ev == categs(ii);
        [Y{ii}, lbins] = countSpikesWithinWindow(sps, t0s(inds), ...
            t1s(inds), tWidth, tShift);
    end
end

function [Y, lbins] = countSpikesWithinWindow(sps, t0s, t1s, ...
    binwidth, binshift)
% 
    if nargin < 5 || isnan(binshift)
        binshift = t1s(1) - t0s(1);
    end
    if nargin < 4 || isnan(binwidth)
        binwidth = t1s(1) - t0s(1);
    end
    lbins = binEdges(0, max(t1s - t0s), binwidth, binshift);
    nbins = numel(lbins);
    Y = nan(numel(t0s), nbins);
    
    for ii = 1:numel(t0s)
        t0 = t0s(ii); t1 = t1s(ii);
        [L, R] = binEdges(t0, t1, binwidth, binshift);
        for jj = 1:numel(L)
            Y(ii,jj) = sum(sps >= L(jj) & sps < R(jj));
        end
    end
end

function [L, R] = binEdges(t0, t1, binwidth, binshift)
    L = t0:binshift:t1-binwidth;
    R = L + binwidth;
end
