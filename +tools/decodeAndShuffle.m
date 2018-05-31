function scs = decodeAndShuffle(X, Y, nshuffles, G)
% for each entry in vs,
%   decode X using Y
%   and decode X using Y shuffled conditional on X
%   if C is provided, shuffle only within each group of C
% 
    if nargin < 4
        G = cell(numel(X), 1);
    end
    nreps = nshuffles; % reps applied only to scsRaw
    scsRaw = nan(numel(X), nreps);
    scsShuf = nan(numel(X), nshuffles);
    nfolds = 5;
    
    parfor ii = 1:numel(X)
        if mod(ii, 10) == 0
            disp(sprintf('Pair %d of %d...', ii, numel(X)));
        end
        
        % load data
        x = X{ii};
        ys = Y{ii};
        gs = G{ii};
        if isempty(gs)
            gs = x == 1;
        end
        ix = ~any(isnan(ys), 2) & ~isnan(x);
        ys = ys(ix,:); x = x(ix); gs = gs(ix);
        if numel(x) == 0
            continue;
        end
        
        % unshuffled
        scsRaw(ii,:) = getDecodingAccuracyOnHoldout(ys, x, nfolds, nreps);
        
        % shuffled
        for jj = 1:nshuffles
%             ix = x == 1;
%             xsp = [x(ix,:); x(~ix,:)];
%             ysp = [randPermByRow(ys(ix,:)); randPermByRow(ys(~ix,:))];
            ysp = shuffleByGroup(ys, gs);
            scsShuf(ii,jj) = getDecodingAccuracyOnHoldout(ysp, ...
                x, nfolds, 1); % only 1 rep
        end
    end
    
    clear scs;
    scs.scsRaw = scsRaw;
    scs.scsShuf = scsShuf;
    scs.scsRawMean = mean(scsRaw, 2);
    scs.scsShufMean = mean(scsShuf, 2);
    scs.scsRawSe = (std(scsRaw,[],2)/sqrt(size(scsRaw,2)));
    scs.scsShufSe = (std(scsShuf,[],2)/sqrt(size(scsShuf,2)));    
    scs.scsDelta = scs.scsRawMean - scs.scsShufMean;
    scs.nreps = nreps;
    scs.nfolds = nfolds;
    scs.nshuffles = nshuffles;
    
end

function ys = shuffleByGroup(y, gs)
    grps = unique(gs);
    ys = nan(size(y));
    for ii = 1:numel(grps)
        ix = gs == grps(ii);
        ys(ix,:) = randPermByRow(y(ix,:));
    end
end

function Yshuf = randPermByRow(Y)
    % randomly shuffle each column of Y1
    [~, idx] = sort(rand(size(Y)),1);
    Yshuf = nan(size(Y));
    for ii = 1:size(Y,2)
        Yshuf(:,ii) = Y(idx(:,ii),ii);
    end
end

function sc = getDecodingAccuracyOnHoldout(Y, X, nfolds, nreps)
    ix = ~any(isnan(Y), 2) & ~isnan(X); Y = Y(ix,:); X = X(ix);
    if numel(X) == 0
        sc = nan;
        return;
    end
    scs = estimate(Y, X, nfolds, nreps);
    sc = scs;
%     sc = nanmean(scs,2); % avg over folds
end

function scs = estimate(Y, X, nfolds, nreps)
    scs = nan(nreps, 1);
    for jj = 1:nreps
%         scs(jj) = rand; continue;
        mdl = fitcdiscr(Y, X, 'CrossVal', 'on', 'KFold', nfolds, ...
            'DiscrimType', 'linear');
        scs(jj) = 1-kfoldLoss(mdl); % accuracy = 1 - err
        continue;
%         ix = crossvalind('Kfold', X, nfolds);
%         for ii = 1:nfolds
%             test = (ix == ii); train = ~test;
%             C = classify(Y(test,:), Y(train,:), X(train), 'linear');
%             scs(jj,ii) = mean(X(test) == C);
%         end
    end
end
