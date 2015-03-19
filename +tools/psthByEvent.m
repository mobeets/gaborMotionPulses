function [Z, lbins, rbins, categs] = psthByEvent(sps, splitEvent, ...
    alignEvent, tL, tR, tWidth, tShift)
%
    inds = ~isnan(alignEvent) & ~isnan(splitEvent);
    alignEvent = alignEvent(inds);
    splitEvent = splitEvent(inds);
    
    t0 = alignEvent - tL;
    t1 = alignEvent + tR;
    [Y, nY, categs] = splitSpikes(sps, t0, t1, splitEvent);
    Z = cell(numel(Y), 1);
    for ii = 1:numel(Y)
        [z, lbins, rbins] = tools.countSpikesWithinWindow(Y{ii}, 0.0, ...
            tR, tWidth, tShift);
        Z{ii} = z./nY(ii);
    end
end

function [Y, nY, categs] = splitSpikes(sps, ts0, ts1, ev)
    categs = sort(unique(ev));
    Y = cell(numel(categs), 1);
    nY = zeros(numel(categs), 1);
    for ii = 1:numel(ts0)
        t0 = ts0(ii); t1 = ts1(ii);        
        ix = sps >= t0 & sps <= t1;
        cind = find(categs == ev(ii));
        Y{cind} = [Y{cind}; sps(ix)-t0];
        nY(cind) = nY(cind) + 1;
    end
    for ii = 1:numel(categs)
        Y{cind} = sort(Y{cind});
    end
end
