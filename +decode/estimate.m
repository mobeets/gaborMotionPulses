function scs = estimate(X, Y, scoreFcn, nfolds, nshuffles)
    scs = nan(nshuffles, nfolds);
    for jj = 1:nshuffles
        ix = crossvalind('Kfold', Y, nfolds);
        for ii = 1:nfolds
            test = (ix == ii); train = ~test;
            [b, dev, stats] = glmfit(X(train,:), Y(train), ...
                'binomial', 'logit');
            Yh = glmval(b, X(test,:), 'logit');
            scs(jj,ii) = scoreFcn(Y(test), round(Yh));
        end
    end
end
