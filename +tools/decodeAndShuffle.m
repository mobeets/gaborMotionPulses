function [scsRaw, scsShuf] = decodeAndShuffle(X, Y, nshuffles, G)
% for each entry in vs,
%   decode X using Y
%   and decode X using Y shuffled conditional on X
%   if C is provided, shuffle only within each group of C
% 
    if nargin < 4
        G = cell(numel(X), 1);
    end
    scsRaw = nan(numel(X), 1);
    scsShuf = nan(numel(X), nshuffles);
    nfolds = 5; nreps = 1;
    
    for ii = 1:numel(X)
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
        scsRaw(ii) = getDecodingAccuracyOnHoldout(ys, x, nfolds, nreps);
        
        % shuffled
        for jj = 1:nshuffles
%             ix = x == 1;
%             xsp = [x(ix,:); x(~ix,:)];
%             ysp = [randPermByRow(ys(ix,:)); randPermByRow(ys(~ix,:))];
            ysp = shuffleByGroup(ys, gs);
            scsShuf(ii,jj) = getDecodingAccuracyOnHoldout(ysp, x, ...
                nfolds, nreps);
        end
    end
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
    sc = nanmean(nanmean(scs,2)); % avg over reps
end

function scs = estimate(Y, X, nfolds, nreps)
    scs = nan(nreps, 1);
    for jj = 1:nreps
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
