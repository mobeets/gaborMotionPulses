function plotCellRF_fit(xy, wf, cmin, sz)
    if nargin < 3 || isnan(cmin)
        cmin = nan;
    end
    if nargin < 4 || isnan(sz)
        sz = 20;
    end
    clrs = plot.getColors(wf, true, cmin);

    hold on;
    axis off;
    axis equal;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        plot(xy(ii,1), xy(ii,2), 'ko', 'MarkerSize', sz, ...
            'MarkerFaceColor', clrs(ii,:));%, 'LineStyle', 'none');
    end
    plot.setPrintSize(gcf, struct('width', 4, 'height', 4));
end
