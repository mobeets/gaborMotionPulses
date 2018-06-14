function fig = dPrimeAndCP(cells, kind)
% kind [str] is either 'dPrime' or 'CP'

    fig = plot.init;
    set(gca, 'LineWidth', 2);
    
    xs = [cells.rf_ecc];
    
    if strcmpi(kind, 'dPrime')    
        ys = [cells.dPrime];
        ylabel('d''');
    elseif strcmpi(kind, 'CP')
        ys = [cells.cp_YresARc];
        ix = ys > 0;
        xs = xs(ix); ys = ys(ix);
        plot([min(xs) max(xs)], [0.5 0.5], 'k--', 'LineWidth', 2);
        ylabel('CP');
    end

    mdl = fitlm(xs, ys);
    h = mdl.plot;
    h(1).Marker = 'o';
    h(1).Color = 'k';
    h(1).MarkerSize = 4.5;
    h(1).MarkerFaceColor = 0.5*ones(3,1);
    h(1).LineWidth = 1.5;
    h(2).LineWidth = 3;
    h(2).Color = [0.8 0.2 0.2];
    set(gca, 'TickDir', 'out');
    legend off;
    xlabel('RF eccentricity');
    xlim([0 1.1]);
    title('');
    plot.setPrintSize(gcf, struct('width', 4, 'height', 3));

end
