function plotKernelSingle(xy, wf, ~, sz, cmin)
    if nargin < 5
        cmin = nan;
    end
    if nargin < 4 || isnan(sz)
        sz = 200;
    end
    clrs = plot.getColors(wf, true, cmin);

    hold on;
%     axis off;
    axis equal;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clrs(ii,:), 'LineStyle', 'none');
    end
end
