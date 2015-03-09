function Y = choiceSelectivity(dt, cellind, stimStartOffset, stimWindowLength, eventName)
    if nargin < 5
        eventName = 'targchosen';
    end
    if nargin < 4
        stimWindowLength = 0.9;
    end
    if nargin < 3
        stimStartOffset = 0.1;
    end
    data = io.loadDataByDate(dt);
    event = data.stim.(eventName);
    neuron = data.neurons{cellind};
    Y = spikesInWindowSplitByEvent(data, event, neuron, ...
        stimStartOffset, stimWindowLength);
    
    bins = linspace(0, max(cell2mat(Y)), 10);
    figure; hold on;
    for ii = 1:numel(Y)
        subplot(numel(Y), 1, ii);
        hist(Y{ii}, bins);
    end
end

function Y = spikesInWindowSplitByEvent(data, event, neuron, t0, t1)
    stim = data.stim;
    inds = stim.goodtrial; % trials without broken fixation
    noNetMotion = sum(sum(stim.pulses, 3), 2) == 0;
    inds = inds & noNetMotion;
    
    motionStartTimes = [stim.timing.motionon] + [stim.timing.plxstart];
    inds = inds & ~isnan(motionStartTimes)';
    ts0 = motionStartTimes(inds) + t0;
    ts1 = ts0 + t1; % const vector of motion length
    
    sps = neuron.spikeTimes;
    Y = countSpikes(sps, ts0, ts1, event(inds));
end

function Y = countSpikes(sps, ts0, ts1, ev)
    categs = sort(unique(ev));    
    Y = cell(numel(categs), 1);
    for ii = 1:numel(ts0)
        t0 = ts0(ii); t1 = ts1(ii);        
        nspikes = sum(sps >= t0 & sps <= t1);
        cind = find(categs == ev(ii));
        Y{cind} = [Y{cind}; nspikes];
    end
end
