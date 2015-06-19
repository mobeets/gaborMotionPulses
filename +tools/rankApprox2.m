function [ms, lbs, ubs, scs] = rankApprox2(fit, data, foldinds, llstr)

    isLinReg = ~strcmp(llstr, 'bern');
    [X, Y, foldinds] = tools.dropTrialsIfYIsNan(data.X, data.Y, foldinds);

    % score functions
    if ~isLinReg
        scoreObj = reg.getScoreObj(isLinReg, 'pctIncorrect');
    else
        scoreObj = reg.getScoreObj(isLinReg, 'rss');
    end
    nullScoreObj = reg.getScoreObj(isLinReg, 'tss');

    % rank-approx functions    
    rankFcn = @(mu, k) reshape(tools.rankKApproximation(...
        reshape(mu, fit.shape(1), fit.shape(2)), k), prod(fit.shape), 1);
    threshMu = @(mu) mu.*(sign(mu)+1)/2; % only keep positive weights

    % fit O.G. and rank-approx models
    ML = reg.getObj_ML(X, Y, struct('foldinds', foldinds));
    scsML = reg.cvFitScores(X, Y, ML, scoreObj);
    Flat = reg.getObj_Flat(X, Y, struct('foldinds', foldinds));
    scsFlat = reg.cvFitScores(X, Y, Flat, scoreObj);
    
    ASD = reg.getObj_ASD(X, Y, data.D, nan, ...
        struct('foldinds', foldinds, 'hyper', fit.hyper));
    fitFcn = ASD.fitFcn;    
    
    % null model
    scsNull = reg.cvFitScores(X, Y, ASD, nullScoreObj);

    % rank-approx models
    rnks = [1 3 fit.shape(2)];
    scs = nan(10, numel(rnks));
    for ii = 1:numel(rnks)
        ASD.fitFcn = @(X, Y, hyper, D) ...
            rankFcn(fitFcn(X,Y,hyper,D), rnks(ii));
        sc = reg.cvFitScores(X, Y, ASD, scoreObj);
        scs(:,ii) = sc;
    end
    
    % summaries
    scsFlatDelta = scs - repmat(scsFlat, 1, size(scs, 2));
    scsNullDelta = scs - repmat(scsNull, 1, size(scs, 2));
    scsDelta = [scsFlatDelta scsNullDelta ...
        scsML-scsFlat scsML-scsNull ...
        scs(:,3)-scs(:,1) scs(:,2)-scs(:,1) scs(:,3)-scsML];

    ms = mean(scsDelta);
    ses = std(scsDelta)/sqrt(size(scsDelta,1));
    lbs = ms - 2*ses;
    ubs = ms + 2*ses;
    
    scs = [scsFlat scsNull scs scsML];
 
end
