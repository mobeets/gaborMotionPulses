function vals = fitSummaries(fitstr)
    vals = struct([]);
    dts = io.getDates('fits');    
    for ii = 1:numel(dts)
        vals = [vals fitSummariesByDate(dts{ii}, fitstr)];
    end
end
