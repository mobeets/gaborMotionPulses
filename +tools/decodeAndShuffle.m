function [scsRaw, scsShuf] = decodeAndShuffle(X, Y, nshuffles)
% for each entry in vs,
%   decode X using Y
%   and decode X using Y shuffled conditional on X
% 
    scsRaw = nan(numel(X), 1);
    scsShuf = nan(numel(X), nshuffles);
    nfolds = 5; nreps = 1;
    
    for ii = 1:numel(X)
        if mod(ii, 10) == 0
            disp(sprintf('Pair %d of %d...', ii, numel(X)));
        end
        x = X{ii};
        ys = Y{ii};
        
        scsRaw(ii) = getDecodingAccuracyOnHoldout(ys, x, nfolds, nreps);
                
        for jj = 1:nshuffles
            ix = x == 1;
            xsp = [x(ix,:); x(~ix,:)];
            ysp = [randPermByRow(ys(ix,:)); randPermByRow(ys(~ix,:))];
            scsShuf(ii,jj) = getDecodingAccuracyOnHoldout(ysp, xsp, ...
                nfolds, nreps);
        end
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
