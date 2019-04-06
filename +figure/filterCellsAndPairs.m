function ix = filterCellsAndPairs(objs, ignoreFlatSpatialRFs, rsqThresh)
    if nargin < 2
        ignoreFlatSpatialRFs = false;
    end
    if nargin < 3
        rsqThresh = -inf;
    end

    ix = true(size(objs));
    ix = ix & ([objs.pctCorrect] >= 0.70);
    ix = ix & ([objs.ntrials] >= 100);
    if isfield(objs, 'dPrime')
        ix = ix & (abs([objs.dPrime]) >= 0.1);
        ix = ix & ~isnan([objs.dPrime]);
    else
        ix = ix & (abs([objs.minAbsDprime]) >= 0.1);
        ix = ix & ~isnan([objs.cell1_dPrime]) & ...
            ~isnan([objs.cell2_dPrime]);
    end
    
    if ignoreFlatSpatialRFs
        if isfield(objs, 'minRfSpatialVar')
            isSpatialVar = log([objs.minRfSpatialVar]) > -15;
        else
            isSpatialVar = log([objs.rfSpatialVariability]) > -15;
        end
        ix = ix & isSpatialVar;
    end
    
    if isfield(objs, 'rsq')
        ix = ix & ([objs.rsq] > rsqThresh);
    else
        ix = ix & ([objs.minRsq] > rsqThresh);
    end

end
