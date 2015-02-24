function [seps, scs] = separabilityVsScore(fitdir)
    dts = {'20130502', '20130514', '20130515', '20130517', '20130611', ...
        '20140213', '20140218', '20140226', '20140303', '20140304', ...
        '20140305', '20140306', '20140307', '20140310'};
    foldind = 1;
    seps = [];
    scs = [];
    for ii = 1:numel(dts)
        dt = dts{ii};
        data = io.loadDataByDate(dt);
        fits = io.loadFitsByDate(dt, fitdir);
        fns = fieldnames(fits);
        for jj = 1:numel(fns)
            fnm = char(fns(jj));
            fit = fits.(fnm);
            [sep, sc] = getVals(fit, data, foldind);
            seps = [seps sep];
            scs = [scs sc];
        end
    end
    figure; plot(seps, scs, 'ok');
    xlabel('separability index'); ylabel('fit score');
end

function [sep, sc] = getVals(fit, data, foldind)
    obj = fit.ASD{end}; % use most recent fit    
    mu = getMu(obj, data.ns, data.nt, foldind);
    sep = getSeparability(mu);
    sc = obj.scores(foldind);
end

function mu = getMu(obj, ns, nt, foldind)
    mu = reshape(obj.mus{foldind}(1:end-1), ns, nt);
end

function val = getSeparability(mu)    
    [~,s,~] = svd(mu);
    s = diag(s);
    val = s(1)/sum(s);
end
