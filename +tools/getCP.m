function [cp, Y, lbins] = getCP(dt, cellind, splitEvent, tL, tR, ...
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
        splitEvent = 'targchosen';
    end
    
    data = io.loadDataByDate(dt);
    stim = data.stim;
    sps = data.neurons{cellind}.spikeTimes;
    splitEvent = stim.(splitEvent);
    alignEvent = [stim.timing.motionon] + [stim.timing.plxstart];
    noNetMotion = sum(sum(stim.pulses, 3), 2) == 0;
    inds = stim.goodtrial & noNetMotion;
    [cp, lbins, Y] = tools.CP(sps, splitEvent, alignEvent', tL, tR, ...
        tWidth, tShift, inds);
end
