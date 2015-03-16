function fig = separabilityVsScore(fitdir)
    dts = io.getDates(fitdir);
    foldind = 1;
    vals = struct([]);
    for ii = 1:numel(dts)
        val = io.summaryByDate(dts{ii}, fitdir, foldind);
        vals = [vals val];
    end
    
    fig = figure; hold on;
    cellTypes = unique({vals.cellType});
    clrs = lines(numel(cellTypes)+1);
    for ii = 1:numel(cellTypes)
        lbl = cellTypes(ii);
        inds = strcmp({vals.cellType}, lbl);
        plot([vals(inds).separability], [vals(inds).score], 'o', ...
            'Color', clrs(ii,:), 'MarkerFaceColor', clrs(ii,:), ...
            'DisplayName', lbl);
    end
    plot([min([vals.separability]) max([vals.separability])], [0 0], ...
        '--', 'Color', 'k', 'HandleVisibility', 'off');
    legend('Location', 'NorthEastOutside');
    xlabel('separability index');
    ylabel('fit score');
end
