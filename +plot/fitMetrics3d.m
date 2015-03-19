function vals = fitMetrics3d(vals, x, y, z)
    figure; hold on;
    xs = [vals.(x)];
    ys = [vals.(y)];
    zs = [vals.(z)];
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    for ii = 1:numel(categs)
        ind = strcmp(categs(ii), {vals.type});
        scatter3(xs(ind), ys(ind), zs(ind), 40, clrs(ii,:), 'filled', ...
            'DisplayName', categs(ii));
    end
    xlabel('score'); ylabel('mu-corr'); zlabel('score-sdev');
    legend(categs, 'Location', 'NorthEastOutside');
end
