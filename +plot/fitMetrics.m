function fitMetrics(vals, x, y, z)
    if nargin < 4
        z = x;
    end
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    figure; hold on;
    
    xs = [vals.(x)];
    ys = [vals.(y)];
    zs = [vals.(z)];
    for ii = 1:numel(categs)
        ind = strcmp({vals.type}, categs(ii));
        scatter3(xs(ind), ys(ind), zs(ind), 40, clrs(ii,:), 'filled', ...
            'DisplayName', categs(ii));
    end
    xlabel(x); ylabel(y); zlabel(z);
    legend(categs, 'Location', 'NorthEastOutside');
%     plot([0 1], [0 1], 'k--');
end
