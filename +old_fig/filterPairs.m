function ix = filterPairs(objs, ignoreFlatSpatialRFs)

    ix = true(size(objs));
    ix = ix & ([objs.pctCorrect] >= 0.70);
    ix = ix & ([objs.ntrials] >= 100);
    
    if isfield(objs, 'minAbsDPrime')
        ix = ix & ([objs.minAbsDprime] >= 0.1);
    else
        ix = ix & ([objs.dPrime] >= 0.1);
    end
    
    isSpatialVar = log([objs.minRfSpatialVar]) > -15;
    if ignoreFlatSpatialRFs
        ix = ix & isSpatialVar;
    end

end
