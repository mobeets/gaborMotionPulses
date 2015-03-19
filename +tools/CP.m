function [cp, ts, Y] = CP(stim, neuron, splitEvent, tL, tR, ...
    tWidth, tShift)
% 
    if nargin < 7
        tShift = nan;
    end
    if nargin < 6
        tWidth = nan;
    end
    if nargin < 5
        tR = 1.35;
    end
    if nargin < 4
        tL = 0.0;
    end    
    if nargin < 3
        splitEventStr = 'targchosen';
        splitEvent = stim.(splitEventStr);
    end
    
    sps = neuron.spikeTimes;
    alignEvent = [stim.timing.motionon] + [stim.timing.plxstart];
    noNetMotion = abs(sum(sum(stim.pulses, 3), 2)) < 10;
    inds = stim.goodtrial & noNetMotion;    
    [Y, ts] = tools.splitAndCountSpikes(sps, splitEvent, ...
        alignEvent', tL, tR, tWidth, tShift, inds);
    % if splitEvent is logical vector, then Y{2} are pref trials
    cp = tools.AUC(Y{2}, Y{1});
end
