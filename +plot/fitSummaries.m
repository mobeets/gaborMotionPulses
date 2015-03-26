function fitSummaries(vals, x, y, z, plotUnity)
    if nargin < 5
        plotUnity = false;
    end
    if nargin < 4 || (numel(z) == 1 && isnan(z))
        z = x;
    end
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    figure; hold on;
    
    xs = [vals.(x)];
    ys = [vals.(y)];
    zs = [vals.(z)];
    for ii = 1:numel(categs)
        if numel(vals) > 1
            ind = strcmp({vals.type}, categs(ii));
        else
            ind = true(numel(xs),1);
        end
        scatter3(xs(ind), ys(ind), zs(ind), 40, clrs(ii,:), 'filled', ...
            'DisplayName', categs(ii));
    end
    lbl = @(x) strrep(x, '_', '-');
    xlabel(lbl(x)); ylabel(lbl(y)); zlabel(lbl(z));
    legend(categs, 'Location', 'NorthEastOutside');
    if plotUnity
        plot([0 1], [0 1], 'k--');
    end
end
