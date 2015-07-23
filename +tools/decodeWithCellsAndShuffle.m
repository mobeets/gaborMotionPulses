function scs = decodeWithCellsAndShuffle(vs, nshuffles)
% for each entry in vs,
%   decode X using Y1 and Y2,
%   and decode X using Y1 and multiple shuffles of Y2
%       where shuffles are done within the same X
% 
    scoreFcn = @(Y, Yh) mean(Y == Yh);
    X = {vs.stim};
    Y1 = {vs.cell1_Y};
    Y2 = {vs.cell2_Y};
    
    scs = cell(numel(X), nshuffles);
    for ii = 1:numel(X)
        if mod(ii,10) == 0
            disp([num2str(ii) ' of ' num2str(numel(X))]);
        end
        x = X{ii};
        y1 = Y1{ii};
        y2 = Y2{ii};
        scs{ii,1} = getMeanCellScore([y1 y2], x, scoreFcn, 10, 10);
        for jj = 2:nshuffles
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
    scs = cell2mat(scs);
end

function sc = getMeanCellScore(X, Y, scoreFcn, nfolds, ntimes)
    ix = all(~isnan(X),2) | isnan(Y); X = X(ix,:); Y = Y(ix);
    if numel(Y) == 0
        sc = nan;
        return;
    end
    sc0 = estimate(X, Y, scoreFcn, nfolds, ntimes);
    sc = nanmean(sc0(:));
end

function scs = estimate(X, Y, scoreFcn, nfolds, ntimes)
    scs = nan(ntimes, nfolds);
    for jj = 1:ntimes
        ix = crossvalind('Kfold', Y, nfolds);
        for ii = 1:nfolds
            test = (ix == ii); train = ~test;
            C = classify(X(test,:), X(train,:), Y(train), 'linear');
            scs(jj,ii) = scoreFcn(Y(test), C);
        end
    end
end
