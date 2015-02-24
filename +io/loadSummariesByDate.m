function [vals, data] = loadSummariesByDate(dt, fitdir, foldind)
    data = io.loadDataByDate(dt);
    fits = io.loadFitsByDate(dt, fitdir);
    fns = fieldnames(fits);
    vals = struct([]);
    for jj = 1:numel(fns)
        fnm = char(fns(jj));
        fit = fits.(fnm);
        curfit = fit.ASD{end}; % use most recent fit
        val = getSummary(curfit, data, fnm, dt, foldind);
        vals = [vals val];
    end
end

function obj = getSummary(fit, data, fnm, dt, foldind)
    obj.name = fnm;
    obj.dt = dt;
    obj.mu = getMu(fit, data.ns, data.nt, foldind);
    obj.mu0 = obj.mu(:,1);
    obj.separability = getSeparability(obj.mu);
    obj.score = fit.scores(foldind);
end

function mu = getMu(obj, ns, nt, foldind)
    mu = reshape(obj.mus{foldind}(1:end-1), ns, nt);
end

function val = getSeparability(mu)    
    [~,s,~] = svd(mu);
    s = diag(s);
    val = s(1)/sum(s);
end
