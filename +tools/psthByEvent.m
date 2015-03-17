function [Z, bins] = psthByEvent(stim, neuron, event, stimEventLength, ntbins, ndiv)
    if nargin < 6
        ndiv = 100;
    end
    if nargin < 5
        ntbins = 10;
    end
    if nargin < 4
        stimEventLength = 1.0;
    end
    inds = stim.goodtrial; % trials without broken fixation
    
    motionStartTimes = [stim.timing.motionon] + [stim.timing.plxstart];
    inds = inds & ~isnan(motionStartTimes)';
    ts0 = motionStartTimes(inds);
    ts1 = ts0 + stimEventLength; % const vector of motion length
    
    sps = neuron.spikeTimes;
    [Y, nY] = splitSpikes(sps, ts0, ts1, event(inds));
    Z = cell(numel(Y), 1);
    for ii = 1:numel(Y)
        [z, bins] = smoothSpikes(Y{ii}, ntbins, ndiv, stimEventLength);
        Z{ii} = z./nY(ii);
    end
end

function [Y, nY] = splitSpikes(sps, ts0, ts1, ev)
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

function [Z, bins] = smoothSpikes(Y, ntbins, ndiv, maxT)
    nbins = ceil(maxT*ndiv)-ntbins;
    Z = zeros(nbins,1);
    bins = zeros(nbins,2);
    g = @(x,y) ([x*ntbins+1 (x+1)*ntbins] + y-2)/ndiv;
    f = @(ii) g(floor(ii/ntbins), mod(ii,ntbins));
    for ii = 1:nbins
        bins(ii,:) = f(ii);
        Z(ii) = sum(Y >= bins(ii,1) & Y <= bins(ii,2));
    end
end
