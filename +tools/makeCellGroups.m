function groups = makeCellGroups(allcells, ks)
% e.g., ks = [2, 3, 4] would give you all pairs, triplets, quads
%
% each cellgroup will have the following fields:
%     dt
%     cellnms
%     stimdir
%     Ys
%     Ys_resAR
%     pctCorrect
%     ntrials
%     minAbsDprime
%     minRfSpatialVar
%

    dts = unique({allcells.dt});
    groups = struct(); % groups: pairs, triplets, quads, etc.
    d = 1;
    for ii = 1:numel(dts)
        cells = allcells(strcmp({allcells.dt}, dts{ii}));
        if numel(cells) == 1
            continue;
        end
        
        % get stimulus direction for decoding
        stim = cells(1);
        X = sign(stim.dirstrength); X(X == 0) = nan;
        
        for kk = 1:numel(ks)
            k = ks(kk);
            cellinds = nchoosek(1:numel(cells), k);
            for jj = 1:size(cellinds,1)
                cinds = cellinds(jj,:);
                ccells = cells(cinds);
                
                groups(d).dt = dts{ii};
                groups(d).cellnms = {ccells.name};
                groups(d).stimdir = X;
                groups(d).ntrials = stim.ntrials;
                groups(d).pctCorrect = stim.pctCorrect;
                groups(d).minAbsDprime = min(abs([ccells.dPrime]));
                groups(d).minRfSpatialVar = min([ccells.rfSpatialVariability]);
                groups(d).targPrefs = [ccells.targPref];
                
                groups(d).Ys = cell2mat({ccells.Y});
                groups(d).Ys_resAR = cell2mat({ccells.YresAR});
                d = d+1;
            end
        end
    end
end
