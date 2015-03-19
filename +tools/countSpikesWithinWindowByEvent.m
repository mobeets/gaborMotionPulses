function [Y, lbins, rbins] = countSpikesWithinWindowByEvent(sps, ev, ...
    t0s, t1s, tWidth, tShift)
% 
    categs = sort(unique(ev));
    Y = cell(numel(categs), 1);
    for ii = 1:numel(categs)
        inds = ev == categs(ii);
        [Y{ii}, lbins, rbins] = tools.countSpikesWithinWindow(sps, ...
            t0s(inds), t1s(inds), tWidth, tShift);
    end
end
