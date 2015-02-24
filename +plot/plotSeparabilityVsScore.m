function vals = plotSeparabilityVsScore(fitdir)
    dts = {'20130502', '20130514', '20130515', '20130517', '20130611', ...
        '20140213', '20140218', '20140226', '20140303', '20140304', ...
        '20140305', '20140306', '20140307', '20140310'};
    foldind = 1;
    vals = struct([]);
    for ii = 1:numel(dts)
        val = io.loadSummariesByDate(dts{ii}, fitdir, foldind);
        vals = [vals val];
    end
    figure; plot([vals.separability], [vals.score], 'ok');
    xlabel('separability index'); ylabel('fit score');
end
