function [scs, scsP, scsA] = decodeWithCells(vs, useAllCells, ...
    predictChoice, scoreFcn, nfolds, nshuffles)
% vs - struct array of all data, from tools.makeFitSummaries()
% scoreFcn - function handle @(Y, Yh) ...;
% 
% 1. score for each monkey
% 2. score for each cell
% 3. score for each cell pair
% 
% dt cell1 cell2 cellScore mnkScore
% 
    if nargin < 6
        nshuffles = 20;
    end
    if nargin < 5
        nfolds = 20;
    end
    if nargin < 4
        scoreFcn = @(Y, Yh) mean(Y == Yh);
    end
    if nargin < 3
        predictChoice = false;
    end
    if nargin < 2
        useAllCells = false;
    end

    dts = unique({vs.dt});
    scs = struct(); % singles
    scsP = struct(); % pairs
    scsA = struct(); % all cells
    c = 1; d = 1; e = 1;
    for ii = 1:numel(dts)
        vc = vs(strcmp({vs.dt}, dts{ii}));
        dec = vc(~[vc.isCell]);
        if numel(dec) == 0
            continue;
        end
        cells = vc([vc.isMT]);
        if numel(cells) == 0
            continue;
        end
        Ym = dec.C; % monkey's choice
        if predictChoice
            Y = Ym; % predict monkey's choice
        else
            Y = sign(dec.dirstrength);
            Y(Y == 0) = nan;
            Y = (Y == 1); % correct choice
        end
        mnkScore = scoreFcn(Y, Ym);

        if useAllCells
            Ys = nan(numel(cells),numel(Y));
            for mm = 1:numel(cells)
                Ys(mm,:) = cells(mm).Y;
            end
            scsA(e).dt = dts{ii};
            scsA(e).mnkScore = mnkScore;
            scsA(e).score = getMeanCellScore(Ys', Y, scoreFcn, ...
                nfolds, nshuffles);
            scsA(e).ncells = numel(cells);
            e = e+1;
            continue;
        end
%         Ynm = 'YhAR';
        Ynm = 'Y';

        for jj = 1:numel(cells)
            Xc1 = cells(jj).(Ynm);
            scs(c).dt = dts{ii};
            scs(c).cell = cells(jj).name;
            scs(c).score = getMeanCellScore(Xc1, Y, scoreFcn, ...
                nfolds, nshuffles);
            scs(c).mnkScore = mnkScore;
            scs(c).dPrime = cells(jj).dPrime;            
            ncA = cells(jj).rf_center;
            c = c+1;
            
            if jj == numel(cells)
                continue;
            end
            for kk = jj+1:numel(cells)
                Xc2 = cells(kk).(Ynm);
                
                scsP(d).dt = dts{ii};
                scsP(d).stim = Y;
                scsP(d).cell1_Y = Xc1;
                scsP(d).cell2_Y = Xc2;
                scsP(d).cell1 = cells(jj).name;
                scsP(d).cell2 = cells(kk).name;
                scsP(d).score = getMeanCellScore([Xc1 Xc2], Y, ...
                    scoreFcn, nfolds, nshuffles);
                scsP(d).mnkScore = mnkScore;
                
                ncB = cells(kk).rf_center';
                scsP(d).rfCorr = corr(cells(jj).w(:), cells(kk).w(:));
                scsP(d).rfDist = sqrt((ncB(1)-ncA(1))^2 + ...
                    (ncB(2) - ncA(2))^2);
                max_dist = sqrt(max(max(asd.sqdist.space(cells(jj).Xxy))));
                scsP(d).rfDist_norm = scsP(d).rfDist/max_dist;
                scsP(d).cell1_dPrime = cells(jj).dPrime;
                scsP(d).cell2_dPrime = cells(kk).dPrime;
                scsP(d).cell1_sep = cells(jj).separability_index;
                scsP(d).cell2_sep = cells(kk).separability_index;
                scsP(d).cell1_targPref = cells(jj).targPref;
                scsP(d).cell2_targPref = cells(kk).targPref;
                scsP(d).sameTarg = scsP(d).cell1_targPref == scsP(d).cell2_targPref;
                scsP(d).sameCorr = scsP(d).rfCorr > 0;
                scsP(d).noiseCorrZer = noiseCorr(cells(jj), cells(kk), ...
                    'Yzer', 30);
                scsP(d).noiseCorrLow = noiseCorr(cells(jj), cells(kk), ...
                    'Ylow', 30);
                scsP(d).noiseCorrAR = noiseCorr(cells(jj), cells(kk), ...
                    'YresAR', 30);

                scsP(d).noiseCorrZerRfCorr = scsP(d).noiseCorrZer * scsP(d).rfCorr;
                scsP(d).noiseCorrLowRfCorr = scsP(d).noiseCorrLow * scsP(d).rfCorr;
                scsP(d).noiseCorrZerRfCorr2 = scsP(d).noiseCorrZer * sign(scsP(d).rfCorr);
                scsP(d).noiseCorrLowRfCorr2 = scsP(d).noiseCorrLow * sign(scsP(d).rfCorr);
                
                scsP(d).dPrimeMax = max([scsP(d).cell1_dPrime, scsP(d).cell2_dPrime]);
                scsP(d).dPrimeMean = mean([scsP(d).cell1_dPrime, scsP(d).cell2_dPrime]);
                d = d+1;
            end
        end
    end
    
    if isempty(fieldnames(scs))
        return;
    end
    
    % backfill with single-cell scores
    for ii = 1:numel(scsP)
        sc1 = scs(strcmp({scs.cell}, scsP(ii).cell1)).score;
        sc2 = scs(strcmp({scs.cell}, scsP(ii).cell2)).score;
        scsP(ii).cell1_score = sc1;
        scsP(ii).cell2_score = sc2;
        scsP(ii).singleScoreMax = max([sc1, sc2]);
        scsP(ii).singleScoreMean = mean([sc1, sc2]);
        scsP(ii).pairImprovement = scsP(ii).score - scsP(ii).singleScoreMax;
        scsP(ii).pairImprovePct = scsP(ii).pairImprovement/scsP(ii).singleScoreMax;
        scsP(ii).pairImprovePctLeft = scsP(ii).pairImprovement/(1-scsP(ii).singleScoreMax);
    end

end

function nc = noiseCorr(c1, c2, Ynm, minTrials)
    Y1 = c1.(Ynm);
    Y2 = c2.(Ynm);
    ix = ~isnan(Y1) & ~isnan(Y2);
    if sum(ix) < minTrials
        nc = nan;
        return;
    end
    Y1 = Y1(ix);
    Y2 = Y2(ix);
    nc = corr(Y1, Y2);
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
%             [b, dev, stats] = glmfit(X(train,:), Y(train), ...
%                 'binomial', 'logit');
%             % use classify( , ,group, 'linear');
%             Yh = glmval(b, X(test,:), 'logit');
%             scs(jj,ii) = scoreFcn(Y(test), round(Yh));
        end
    end
end

