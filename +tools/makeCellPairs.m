function pairs = makeCellPairs(allcells)

    dts = unique({allcells.dt});
    pairs = struct(); % pairs
    d = 1;
    for ii = 1:numel(dts)
        cells = allcells(strcmp({allcells.dt}, dts{ii}));
        
        % get stimulus direction for decoding
        %   and also bin dirstrength
        stim = cells(1);
        X = sign(stim.dirstrength); X(X == 0) = nan;
        dirs = stim.dirstrength;
        if numel(unique(dirs)) > 10
            dirsbn = discretize(dirs, linspace(min(dirs), max(dirs), 10));
        else
            dirsbn = dirs;
        end

        for jj = 1:(numel(cells)-1)
            for kk = jj+1:numel(cells)
                pairs(d).dt = dts{ii};
                pairs(d).cell1 = cells(jj).name;
                pairs(d).cell2 = cells(kk).name;
                pairs(d).Ys = [cells(jj).Y cells(kk).Y];
                pairs(d).stimdir = X;
                pairs(d).dirstrength = dirs;
                pairs(d).dirstrength_binned = dirsbn;                
                
                pairs(d).rfCorr = corr(cells(jj).w(:), cells(kk).w(:));
                rf1 = cells(jj).wsep.spatial_RF;
                rf2 = cells(kk).wsep.spatial_RF;
                if sum(cells(jj).wsep.temporal_RF < 0) > 3
                    rf1 = -rf1;
                end
                if sum(cells(kk).wsep.temporal_RF < 0) > 3
                    rf2 = -rf2;
                end
                pairs(d).rfCorr_spatial = corr(rf1, rf2);
                pairs(d).sameCorr_spatial = pairs(d).rfCorr_spatial > 0;
                
                pairs(d).sameCorr = pairs(d).rfCorr > 0;
                pairs(d).rfDist = norm(cells(jj).rf_center - cells(kk).rf_center);
                pairs(d).cell1_dPrime = cells(jj).dPrime;
                pairs(d).cell2_dPrime = cells(kk).dPrime;
                pairs(d).cell1_targPref = cells(jj).targPref;
                pairs(d).cell2_targPref = cells(kk).targPref;
                pairs(d).sameTarg = pairs(d).cell1_targPref == pairs(d).cell2_targPref;
                pairs(d).noiseCorrAR = tools.noiseCorr(cells(jj).YresAR, ...
                    cells(kk).YresAR, 30);
                
                % for later filtering
                pairs(d).ntrials = stim.ntrials;
                pairs(d).pctCorrect = stim.pctCorrect;
                pairs(d).minAbsDprime = min(abs(pairs(d).cell1_dPrime), ...
                    abs(pairs(d).cell2_dPrime));
                pairs(d).minRfSpatialVar = min(cells(jj).rfSpatialVariability, ...
                    cells(kk).rfSpatialVariability);
                d = d+1;
            end
        end
    end
end
