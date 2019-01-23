function ix = filterCellsAndPairs(objs, ignoreFlatSpatialRFs)

    ix = true(size(objs));
    ix = ix & ([objs.pctCorrect] >= 0.70);
    ix = ix & ([objs.ntrials] >= 100);
    if isfield(objs, 'dPrime')
        ix = ix & (abs([objs.dPrime]) >= 0.1);
    else
        ix = ix & (abs([objs.minAbsDprime]) >= 0.1);
    end
    
    if ignoreFlatSpatialRFs
        if isfield(objs, 'minRfSpatialVar')
            isSpatialVar = log([objs.minRfSpatialVar]) > -15;
        else
            isSpatialVar = log([objs.rfSpatialVariability]) > -15;
        end
        ix = ix & isSpatialVar;
    end

end
