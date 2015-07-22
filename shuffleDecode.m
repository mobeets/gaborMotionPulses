function scs = shuffleDecode(vs, npermutes)
    scoreFcn = @(Y, Yh) mean(Y == Yh);
    X = {vs.stim};
    Y1 = {vs.cell1_Y};
    Y2 = {vs.cell2_Y};
    
    scs = cell(numel(X), npermutes);
    for ii = 1:numel(X)
        if mod(ii,10) == 0
            disp([num2str(ii) ' of ' num2str(numel(X))]);
        end
        x = X{ii};
        y1 = Y1{ii};
        y2 = Y2{ii};
        scs{ii,1} = getMeanCellScore([y1 y2], x, scoreFcn, 10, 10);
        for jj = 2:npermutes
            ix0 = x == 1;
            Xc = [x(ix0); x(~ix0)];
            Y1c = [y1(ix0); y1(~ix0)];
            
            ix1 = randperm(sum(ix0));
            ix2 = randperm(sum(~ix0));
            Y2a = y2(ix0);
            Y2a = Y2a(ix1);
            Y2b = y2(~ix0);
            Y2b = Y2b(ix2);
            Y2c = [Y2a; Y2b];
            
            Yc = [Y1c Y2c];
            scs{ii,jj} = getMeanCellScore(Yc, Xc, scoreFcn, 10, 10);
        end
    end
    
end

function sc = getMeanCellScore(X, Y, scoreFcn, nfolds, nshuffles)
    ix = all(~isnan(X),2) | isnan(Y); X = X(ix,:); Y = Y(ix);
    if numel(Y) == 0
        sc = nan;
        return;
    end
    sc0 = estimate(X, Y, scoreFcn, nfolds, nshuffles);
    sc = nanmean(sc0(:));
end

function scs = estimate(X, Y, scoreFcn, nfolds, nshuffles)
    scs = nan(nshuffles, nfolds);
    for jj = 1:nshuffles
        ix = crossvalind('Kfold', Y, nfolds);
        for ii = 1:nfolds
            test = (ix == ii); train = ~test;
            C = classify(X(test,:), X(train,:), Y(train), 'linear');
            scs(jj,ii) = scoreFcn(Y(test), C);
        end
    end
end
