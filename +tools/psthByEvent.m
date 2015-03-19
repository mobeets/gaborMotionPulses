function [Z, lbins, rbins, categs] = psthByEvent(stim, neuron, ...
    splitEvent, windowLength, binwidth, binshift)
% 
    if nargin < 6
        binshift = 0.01;
    end
    if nargin < 5
        binwidth = 0.1;
    end
    if nargin < 4
        windowLength = 1.35;
    end
    sps = neuron.spikeTimes;
    inds = stim.goodtrial; % trials without broken fixation
    inds0 = false(numel(stim.goodtrial),1);
    inds0(neuron.trialIndex) = true;
    inds = inds & inds0;
%     inds = inds & abs(sum(sum(stim.pulses,3),2)) < 10;
    
    motionStartTimes = [stim.timing.motionon] + [stim.timing.plxstart];
    inds = inds & ~isnan(motionStartTimes)';
    ts0 = motionStartTimes(inds);
    ts1 = ts0 + windowLength; % const vector of motion length
    
    [Y, nY, categs] = splitSpikes(sps, ts0, ts1, splitEvent(inds));
    Z = cell(numel(Y), 1);
    for ii = 1:numel(Y)
        [z, lbins, rbins] = tools.countSpikesWithinWindow(Y{ii}, 0.0, ...
            windowLength, binwidth, binshift);
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
