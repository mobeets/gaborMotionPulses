function ix = filterPairs(pairs)

    ix = true(size(pairs));
    ix = ix & ([pairs.pctCorrect] >= 0.70);
    ix = ix & ([pairs.ntrials] >= 100);
    ix = ix & ([pairs.minAbsDprime] >= 0.1);
    
    isSpatialVar = log([pairs.minRfSpatialVar]) > -15;
    ix = ix & isSpatialVar;
    disp(['Keeping ' num2str(sum(ix)) ' cell pairs, discarding ' ...
        num2str(sum(~ix))]);

end
