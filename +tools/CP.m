function [cp, ts, Y] = CP(sps, splitEvent, alignEvent, tL, tR, ...
    tWidth, tShift)
% 
    inds = ~isnan(alignEvent) & ~isnan(splitEvent);
    alignEvent = alignEvent(inds);
    splitEvent = splitEvent(inds);
    [Y, ts] = tools.splitAndCountSpikes(sps, splitEvent, ...
        alignEvent, tL, tR, tWidth, tShift);
    % if splitEvent is logical vector, then Y{2} are pref trials
    cp = tools.AUC(Y{2}, Y{1});
end
