function [scs0, scsSh] = decodeWithCellsAndShuffle(X, Y, nshuffles)
% for each entry in vs,
%   decode X using Y
%   and decode X using Y shuffled conditional on X
% 
    scoreFcn = @(Y, Yh) mean(Y == Yh);
%     X = {vs.stim};
%     Y = {vs.Ys};
    
    scs0 = nan(numel(X), 1);
    scsSh = nan(numel(X), nshuffles);
    for ii = 1:numel(X)
        if mod(ii,10) == 0
            disp([num2str(ii) ' of ' num2str(numel(X))]);
        end
        x = X{ii};
        ys = Y{ii}';
        ix = x == 1;
        x = [x(ix); x(~ix)];
        ys = [ys(ix,:); ys(~ix,:)];
        ix = x == 1;
        scs0(ii) = getMeanCellScore(ys, x, scoreFcn, 10, 10);
        for jj = 1:nshuffles
            ysp1 = randPermByRow(ys(ix,:));
            ysp2 = randPermByRow(ys(~ix,:));
            ysp = [ysp1; ysp2];
            scsSh(ii,jj) = getMeanCellScore(ysp, x, scoreFcn, 10, 10);
        end
    end
end

function Y2 = randPermByRow(Y1)
    [~, idx] = sort(rand(size(Y1)),1);
    Y2 = cell2mat(arrayfun(@(i) Y1(idx(:,i),i), 1:size(Y1,2), 'uni', 0));
end

function sc = getMeanCellScore(X, Y, scoreFcn, nfolds, ntimes)
    ix = all(~isnan(X),2) | isnan(Y); X = X(ix,:); Y = Y(ix);
    if numel(Y) == 0
        sc = nan;
        return;
    end
    [sc0, lgps] = estimate(X, Y, scoreFcn, nfolds, ntimes);
    sc = nanmean(sc0(:));
end

function [scs, lgps] = estimate(X, Y, scoreFcn, nfolds, ntimes)
    scs = nan(ntimes, nfolds);
    lgps = cell(ntimes, nfolds);
    for jj = 1:ntimes
        ix = crossvalind('Kfold', Y, nfolds);
        for ii = 1:nfolds
            test = (ix == ii); train = ~test;
            [C,~,~,lgp] = classify(X(test,:), X(train,:), Y(train), 'linear');
            lgps{jj,ii} = lgp;
            scs(jj,ii) = scoreFcn(Y(test), C);
        end
    end
end
