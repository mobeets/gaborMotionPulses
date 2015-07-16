function v = autoRegressModelSpikes(v, nlags, nfolds, nshuffles)
    if nargin < 3
        nfolds = 10;
    end
    if nargin < 4
        nshuffles = 10;
    end
    scoreFcn = @tools.rsq;
    
    Y = v.Y;
    ix = ~isnan(Y);
    Y0 = Y(ix);
    Yh = v.Yh(ix);
    sc0 = scoreFcn(Yh, Y0);

    [X, Y] = makeLagMats(Yh, Y0, nlags);    
    mdl = fitlm(X, Y);
    Yh1 = mdl.predict(X);
    Yh1 = [Yh(1:nlags); Yh1];
%     scores = estimate(X, Y, scoreFcn, nfolds, nshuffles);
%     sc1 = nanmean(scores(:));
    sc1 = scoreFcn(Yh1, Y0);
%     disp(num2str([sc0 sc1]));
    
    v.YhAR = nan(size(v.Y));
    v.YhAR(ix) = Yh1;
    v.score_AR = sc1;
    v.AR_nlags = nlags;

end

function [X, Y] = makeLagMats(Yh, Y0, nlags)
    nlags = nlags+1;
    X = nan(numel(Y0)-nlags+1, nlags+1);
    X(:,end) = Yh(nlags:end);
    for jj = 1:nlags
        X(:,jj) = Y0(nlags-jj+1:end-jj+1);
    end
    Y = X(:,1);
    X = X(:,2:end);
end

function scs = estimate(X, Y, scoreFcn, nfolds, nshuffles)
    scs = nan(nshuffles, nfolds);
    for jj = 1:nshuffles
        ix = crossvalind('Kfold', Y, nfolds);
        for ii = 1:nfolds
            test = (ix == ii); train = ~test;
            mdl = fitlm(X(train,:), Y(train));
            Yh = mdl.predict(X(test,:));
            scs(jj,ii) = scoreFcn(Yh, Y(test));
        end
    end    
end
