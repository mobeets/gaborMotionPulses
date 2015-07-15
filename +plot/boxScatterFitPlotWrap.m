function boxScatterFitPlotWrap(vs, xlbl, ylbl, sB, sF, nb)
    if nargin < 4
        sB = true;
    end
    if nargin < 5
        sF = true;
    end
    if nargin < 6
        nb = 12;
    end
    figure; hold on;
    set(gcf, 'color', 'w');
    xs = [vs.(xlbl)];
    ys = [vs.(ylbl)];
    ix = ~isnan(xs) & ~isnan(ys);
    xs = xs(ix);
    ys = ys(ix);
    plot.boxScatterFitPlot(xs', ys', sB, sF, nb);
    xlabel(strrep(xlbl, '_', '-'));
    ylabel(strrep(ylbl, '_', '-'));
end
