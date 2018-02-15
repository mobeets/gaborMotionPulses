function pairs = makeCellPairs(allcells)

    dts = unique({allcells.dt});
    pairs = struct(); % pairs
    d = 1;
    for ii = 1:numel(dts)
        cells = allcells(strcmp({allcells.dt}, dts{ii}));
        
        % get stimulus direction for decoding
        dec = cells(1);
        Y = sign(dec.dirstrength);
        Y(Y == 0) = nan;
%         Y = (Y == 1); % correct choice

        for jj = 1:(numel(cells)-1)
            for kk = jj+1:numel(cells)
                pairs(d).dt = dts{ii};
                pairs(d).cell1 = cells(jj).name;
                pairs(d).cell2 = cells(kk).name;
                pairs(d).Ys = [cells(jj).Y cells(kk).Y];
                pairs(d).stim = Y;
                
                pairs(d).rfCorr = corr(cells(jj).w(:), cells(kk).w(:));
                pairs(d).sameCorr = pairs(d).rfCorr > 0;
                pairs(d).rfDist = norm(cells(jj).rf_center - cells(kk).rf_center);
                pairs(d).cell1_dPrime = cells(jj).dPrime;
                pairs(d).cell2_dPrime = cells(kk).dPrime;
                pairs(d).cell1_sep = cells(jj).separability_index;
                pairs(d).cell2_sep = cells(kk).separability_index;
                pairs(d).cell1_targPref = cells(jj).targPref;
                pairs(d).cell2_targPref = cells(kk).targPref;
                pairs(d).sameTarg = pairs(d).cell1_targPref == pairs(d).cell2_targPref;
                pairs(d).noiseCorrAR = tools.noiseCorr(cells(jj).YresAR, ...
                    cells(kk).YresAR, 30);
                d = d+1;
            end
        end
    end
end
