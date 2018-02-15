function [scs0, scsSh] = decodeWithCellsAndShuffle(X, Y, nshuffles, doPlot)
% for each entry in vs,
%   decode X using Y
%   and decode X using Y shuffled conditional on X
% 
    if nargin < 4
        doPlot = false;
    end
    scoreFcn = @(Y, Yh) mean(Y == Yh);
%     X = {vs.stim};
%     Y = {vs.Ys};
    
    scs0 = nan(numel(X), 1);
    scsSh = nan(numel(X), nshuffles);
    L0 = nan(numel(X), 2);
    K0 = nan(numel(X), 1);
    Ls = nan(numel(X), nshuffles, 2);
    Ks = nan(numel(X), nshuffles);
    
    nr = 1; nc = nshuffles+1;
    if nc > 5
        nr = floor(sqrt(nshuffles));
        nc = ceil(nshuffles/nr);
    end
    
    for ii = 1:numel(X)
        if doPlot
            figure(2*ii); clf; hold on;
        end
        if mod(ii,10) == 0
            disp([num2str(ii) ' of ' num2str(numel(X))]);
        end
        x = X{ii};
        ys = Y{ii};
        ix = x == 1;
        x = [x(ix); x(~ix)];
        ys = [ys(ix,:); ys(~ix,:)];
        ix = x == 1;
        [scs0c, L0c, K0c] = getMeanCellScore(ys, x, ...
            scoreFcn, 10, 10);
        
        scs0(ii) = scs0c;
        L0(ii,:) = L0c;
        K0(ii) = K0c;
        if doPlot
            subplot(nr,nc,1); hold on;
            scatter(ys(ix,1), ys(ix,2), 'b');
            scatter(ys(~ix,1), ys(~ix,2), 'r');
        end
        for jj = 1:nshuffles
            ysp1 = randPermByRow(ys(ix,:));
            ysp2 = randPermByRow(ys(~ix,:));
            ysp = [ysp1; ysp2];
            [scsSh(ii,jj), Ls(ii,jj,:), Ks(ii,jj)] = getMeanCellScore(...
                ysp, x, scoreFcn, 10, 10);
            if doPlot
                subplot(nr,nc,jj+1); hold on;
                scatter(ysp(ix,1), ysp(ix,2), 'b');
                scatter(ysp(~ix,1), ysp(~ix,2), 'r');
            end
        end
        
        if doPlot
            L = L0(ii,:)';
            K = K0(ii);
            x0 = ys*L + K;
            figure(2*ii-1); clf;
            subplot(nr,nc,1); hold on;
            minx = min(x0(:));
            maxx = max(x0(:));
            bins = linspace(minx, maxx, 20);
            bs = histc(x0(ix), bins);
            stairs(bins, bs, 'b');
            bs = histc(x0(~ix), bins);
            stairs(bins, bs, 'r');
            title(scs0(ii));

            for jj = 1:nshuffles
                L2 = squeeze(Ls(ii,jj,:));
                K2 = Ks(ii,jj);
                x1 = ys*L2 + K2;
                subplot(nr,nc,jj+1); hold on;
                bs = histc(x1(ix), bins);
                stairs(bins, bs, 'b');
                bs = histc(x1(~ix), bins);
                stairs(bins, bs, 'r');
                title(scsSh(ii,jj));
            end
        end
    end
end

function Y2 = randPermByRow(Y1)
    [~, idx] = sort(rand(size(Y1)),1);
    Y2 = cell2mat(arrayfun(@(i) Y1(idx(:,i),i), 1:size(Y1,2), 'uni', 0));
end

function [sc, L, K] = getMeanCellScore(X, Y, scoreFcn, nfolds, ntimes)
    ix = all(~isnan(X),2) | isnan(Y); X = X(ix,:); Y = Y(ix);
    if numel(Y) == 0
        sc = nan; L = nan; K = nan;
        return;
    end
    sc0 = estimate(X, Y, scoreFcn, nfolds, ntimes);    
    sc = nanmean(sc0(:));
    [C,~,~,~,c] = classify(X, X, Y, 'linear');
    L = c(2).linear; K = c(2).const;
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
