function [Z, lbins, rbins, categs] = getPsthByEvent(stim, neuron, ...
    splitEvent, windowLength, tWidth, tShift)
% 
    if nargin < 6
        tShift = 0.01;
    end
    if nargin < 5
        tWidth = 0.1;
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
    motionStartTimes = motionStartTimes(inds)';
    splitEvent = splitEvent(inds);
    
    [Z, lbins, rbins, categs] = tools.psthByEvent(sps, splitEvent, ...
        motionStartTimes, 0.0, windowLength, tWidth, tShift);
end
