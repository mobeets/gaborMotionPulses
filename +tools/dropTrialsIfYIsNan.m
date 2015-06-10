function [X, Y, foldinds, evalinds] = dropTrialsIfYIsNan(X, Y, ...
    foldinds, evalinds)
    inds = isnan(Y);
    X = X(~inds,:,:);
    Y = Y(~inds,:);
    
    % remove foldinds relative to longer inds list
    if numel(inds) > numel(foldinds)
        f_inds = nan(numel(inds),1);
        f_inds(evalinds) = foldinds;
        f_inds(inds) = NaN;
        foldinds = f_inds(~isnan(f_inds));
    else
        foldinds = foldinds(~inds);
    end
    
    if nargin > 3
        evalinds = evalinds(~inds);
    end
end
