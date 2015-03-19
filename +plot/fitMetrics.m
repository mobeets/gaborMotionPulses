function fitMetrics(vals, x, y)
    figure; hold on;
    xlabel(x); ylabel(y);
    
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    xs = [vals.(x)];
    ys = [vals.(y)];    
    for ii = 1:numel(categs)
        inds = strcmp({vals.type}, categs(ii));
        scatter(xs(inds), ys(inds), 50, ...
            'DisplayName', categs(ii), 'MarkerFaceColor', clrs(ii,:));
    end
    legend(categs, 'Location', 'NorthEastOutside');
%     plot([0 1], [0 1], 'k--');
end
