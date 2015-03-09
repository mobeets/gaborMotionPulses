function fig = separabilityVsScore(fitdir)
    dts = io.getDates(fitdir);
    foldind = 1;
    vals = struct([]);
    for ii = 1:numel(dts)
        val = io.summaryByDate(dts{ii}, fitdir, foldind);
        vals = [vals val];
    end
    fig = figure; plot([vals.separability], [vals.score], 'ok');
    xlabel('separability index'); ylabel('fit score');
end
