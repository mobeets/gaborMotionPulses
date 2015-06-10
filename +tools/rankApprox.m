function [ms, lbs, ubs, scs] = rankApprox(fit, data, foldinds, llstr)

    [MAPScoreFcn, scoreObj] = reg.scoreFcns('rsq', llstr);
    if strcmp(llstr, 'bern')
        nullScoreFcn = scoreObj.pctIncorrectFcn;
        rankScoreFcn = scoreObj.pctIncorrectFcn;
    else
        nullScoreFcn = scoreObj.tssFcn;
        rankScoreFcn = scoreObj.rssFcn;
    end
    MLScoreFcn = rankScoreFcn;
    
    MAPFcn = @(hyper) asd.fitHandle(fit.hyper, data.D, llstr);
    MLFcn = @(~) ml.fitHandle(llstr);
    
    rankFcn = @(mu, k) reshape(tools.rankKApproximation(...
        reshape(mu, fit.shape(1), fit.shape(2)), k), prod(fit.shape), 1);
    
    
    threshMu = @(mu) mu.*(sign(mu)+1)/2; % only keep positive weights
    rankFcns = {@(mu) ones(size(mu,1)-1, 1), ...
        @(mu) threshMu(mu(1:end-1)), ... % n.b. no offset
        @(mu) rankFcn(mu, 1), ...
        @(mu) rankFcn(mu, 3), ...
        @(mu) rankFcn(mu, 7)};
    scoreFcns = [{nullScoreFcn}, ...
        repmat({rankScoreFcn}, 1, numel(rankFcns)-1)];
    
    [X, Y, foldinds] = tools.dropTrialsIfYIsNan(data.X, data.Y, foldinds);
    trials = reg.trainAndTestKFolds(X, Y, nan, foldinds);    

    if size(fit.hyper,2) > 1
        fit.hyper = fit.hyper';
    end

    [~,~,mus0] = reg.cvScoreGrid(trials, MAPFcn, MAPScoreFcn, fit.hyper');
    [scsML,~,~] = reg.cvScoreGrid(trials, MLFcn, MLScoreFcn, nan);
    [scs, ~] = cvRankApprox(trials, rankFcns, scoreFcns, mus0, fit.hyper);
        
    nullSc = scs(:,1);
    scs = scs(:,2:end);
    scsDelta = scs - repmat(nullSc, 1, size(scs, 2));
    scsDelta = [scs(:,2)-scs(:,3) scs(:,2)-scs(:,4) scs(:,1)-scs(:,4) ...
        scsDelta scs(:,4)-scsML' scsML'-nullSc];
    
    ms = mean(scsDelta);
    ses = std(scsDelta)/sqrt(size(scsDelta,1));
    lbs = ms - 2*ses;
    ubs = ms + 2*ses;
    
    scs = [nullSc scsML' scs];
    
end

function [scores, mus] = cvRankApprox(trials, fitFcns, scoreFcns, mus0, hyper)
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
            scoreFcn = scoreFcns{jj};
            wf = fitFcn(wf0);
            mus{ii,jj} = [wf; b];
            scores(ii, jj) = scoreFcn(ctrials, mus{ii,jj}, hyper);
        end
    end
end
