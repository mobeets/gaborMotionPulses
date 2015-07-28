function [nc, ps] = noiseCorr(Y1, Y2, minTrials, C)
    ix = ~isnan(Y1) & ~isnan(Y2);
    if nargin < 4
        C = ones(numel(Y1),1);
    end
    conds = unique(C);
    nc = nan(numel(conds),1);
    ps = nan(numel(conds),1);
    for ii = 1:numel(conds)
       idx = ix & (C == conds(ii));
       if sum(idx) < minTrials
           nc(ii) = nan;
           ps(ii) = nan;
       else
           [nc(ii), ps(ii)] = corr(Y1(idx), Y2(idx));
       end
    end
end
