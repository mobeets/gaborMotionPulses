function [ms, lbs, ubs, scs] = rankApprox2(fit, data, foldinds, llstr)

    isLinReg = ~strcmp(llstr, 'bern');
    [X, Y, foldinds] = tools.dropTrialsIfYIsNan(data.X, data.Y, foldinds);
    trials = tools.trainAndTestKFolds(X, Y, nan, foldinds);
%     if size(fit.hyper,2) > 1
%         fit.hyper = fit.hyper';
%     end


    % score functions
    if ~isLinReg
        scoreObj = reg.getScoreObj(isLinReg, 'pctIncorrect');
    else
        scoreObj = reg.getScoreObj(isLinReg, 'rss');
    end

    % rank-approx functions    
    rankFcn = @(mu, k) reshape(tools.rankKApproximation(...
        reshape(mu, fit.shape(1), fit.shape(2)), k), prod(fit.shape), 1);
    
    threshMu = @(mu) mu.*(sign(mu)+1)/2; % only keep positive weights
    rankFcns = {@(mu) threshMu(mu(1:end-1)), ... % n.b. no offset
        @(mu) rankFcn(mu, 1), ...
        @(mu) rankFcn(mu, 3), ...
        @(mu) rankFcn(mu, 7)};
    scoreFcns = repmat({scoreObj}, 1, numel(rankFcns));

    % fit O.G. and rank-approx models
    ASD = reg.getObj_ASD(X, Y, data.D, nan, ...
        struct('foldinds', foldinds, 'hyper', fit.hyper));
    [~,mus0,~,~] = reg.cvFitScores(X, Y, ASD, scoreObj);
    ML = reg.getObj_ML(X, Y, struct('foldinds', foldinds));
    [scsML,~,~,~] = reg.cvFitScores(X, Y, ML, scoreObj);
    Flat = reg.getObj_Flat(X, Y, struct('foldinds', foldinds));
    [scsFlat,~,~,~] = reg.cvFitScores(X, Y, Flat, scoreObj);
    
    fitFcn = ASD.fitFcn;
    scs = [];
    ASD.fitFcn = @(X, Y, hyper, D) ...
        threshMu(fitFcn(X,Y,hyper,D));
    [sc, ~] = reg.cvFitScores(X, Y, ASD, scoreObj);
    scs = [scs sc];
    for ii = [1 3 7]
        ASD.fitFcn = @(X, Y, hyper, D) ...
            rankFcn(fitFcn(X,Y,hyper,D), ii);
        [sc, ~] = reg.cvFitScores(X, Y, ASD, scoreObj);
        scs = [scs sc];
    end
    
    % summaries
    scsDelta = scs - repmat(scsFlat, 1, size(scs, 2));
    scsDelta = [scs(:,2)-scs(:,3) scs(:,2)-scs(:,4) scs(:,1)-scs(:,4) ...
        scsDelta scs(:,4)-scsML scsML-scsFlat];
    
    ms = mean(scsDelta);
    ses = std(scsDelta)/sqrt(size(scsDelta,1));
    lbs = ms - 2*ses;
    ubs = ms + 2*ses;
    
    scs = [scsFlat scsML scs];
 
end
