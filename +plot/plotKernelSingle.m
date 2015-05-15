function plotKernelSingle(xy, wf, vmax)
    if nargin < 3
        vmax = max(abs(wf(:)));
    end
    sz = 220;
    clrFcn = plot.colorScheme();
    wf = wf/vmax;
    hold on;
    axis off; axis square;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        clr = clrFcn(wf(ii));
        plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clr, 'LineStyle', 'none');
    end
end
