function [Y, lbins, rbins] = countSpikesWithinWindow(sps, t0s, t1s, ...
    binwidth, binshift)
% 
    if nargin < 5 || isnan(binshift)
        binshift = t1s(1) - t0s(1);
    end
    if nargin < 4 || isnan(binwidth)
        binwidth = t1s(1) - t0s(1);
    end
    [lbins, rbins] = binEdges(0, max(t1s - t0s), binwidth, binshift);
    nbins = numel(lbins);
    Y = nan(numel(t0s), nbins);
    
    for ii = 1:numel(t0s)
        t0 = t0s(ii); t1 = t1s(ii);
        [L, R] = binEdges(t0, t1, binwidth, binshift);
        if numel(L) > numel(lbins)
            % n.b. precision errors can lead to lbins being too short
            lbins = L - t0;
            rbins = R - t0;
        end
        for jj = 1:numel(L)
            Y(ii,jj) = sum(sps >= L(jj) & sps < R(jj));
        end
    end
end

function [L, R] = binEdges(t0, t1, binwidth, binshift)
    L = t0:binshift:t1-binwidth;
    R = L + binwidth;
end
