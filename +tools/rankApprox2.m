function [ms, lbs, ubs, scs] = rankApprox2(fit, data, foldinds, llstr)

    isLinReg = ~strcmp(llstr, 'bern');
    [X, Y, foldinds] = tools.dropTrialsIfYIsNan(data.X, data.Y, foldinds);
    trials = tools.trainAndTestKFolds(X, Y, nan, foldinds);
%     if size(fit.hyper,2) > 1
%         fit.hyper = fit.hyper';
%     end


    % score functions
    scoreObj = reg.getScoreObj(isLinReg, 'rsq');    
    if ~isLinReg
%         nullScoreFcn = reg.getScoreObj(isLinReg, 'pctIncorrect');
        rankScoreFcn = reg.getScoreObj(isLinReg, 'pctIncorrect');
    else
%         nullScoreFcn = reg.getScoreObj(isLinReg, 'tss');
        rankScoreFcn = reg.getScoreObj(isLinReg, 'rsq');
    end
    scoreObjML = rankScoreFcn;
    scoreObjFlat = rankScoreFcn;

    % rank-approx functions    
    rankFcn = @(mu, k) reshape(tools.rankKApproximation(...
        reshape(mu, fit.shape(1), fit.shape(2)), k), prod(fit.shape), 1);
    
    threshMu = @(mu) mu.*(sign(mu)+1)/2; % only keep positive weights
    rankFcns = {@(mu) threshMu(mu(1:end-1)), ... % n.b. no offset
        @(mu) rankFcn(mu, 1), ...
        @(mu) rankFcn(mu, 3), ...
        @(mu) rankFcn(mu, 7)};
    scoreFcns = repmat({rankScoreFcn}, 1, numel(rankFcns));

    % fit O.G. and rank-approx models
    ASD = reg.getObj_ASD(X, Y, data.D, nan, ...
        struct('foldinds', foldinds, 'hyper', fit.hyper));
    [~,mus0,~,~] = reg.cvFitScores(X, Y, ASD, scoreObj);
    ML = reg.getObj_ML(X, Y, struct('foldinds', foldinds));
    [scsML,~,~,~] = reg.cvFitScores(X, Y, ML, scoreObjML);
    Flat = reg.getObj_Flat(X, Y, struct('foldinds', foldinds));
    [scsFlat,~,~,~] = reg.cvFitScores(X, Y, Flat, scoreObjFlat);
    [scs, ~] = cvRankApprox(trials, rankFcns, scoreFcns, mus0);
    
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

function [scores, mus] = cvRankApprox(trials, fitFcns, sObjs, mus0)
    nfolds = numel(trials);
    nfcns = numel(fitFcns);
    scores = nan(nfolds, nfcns);
    mus = cell(nfolds, nfcns);
    for ii = 1:nfolds
        ctrials = trials(ii);
        mu0 = mus0{ii};
        wf0 = mu0(1:end-1);
        b = mu0(end);
        for jj = 1:nfcns
            fitFcn = fitFcns{jj};
            so = sObjs{jj};
            wf = fitFcn(wf0);
            mus{ii,jj} = [wf; b];
            scores(ii, jj) = so.scoreFcn(ctrials.x_train, ...
                ctrials.y_train, mus{ii,jj}, ctrials);
        end
    end
end
