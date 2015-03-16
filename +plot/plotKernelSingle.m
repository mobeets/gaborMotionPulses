function plotKernelSingle(xy, wf, vmax)
    if nargin < 3
        vmax = max(abs(wf(:)));
    end
    sz = 50;
    clrFcn = plot.colorScheme();
    wf = wf/vmax;
    hold on;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        clr = clrFcn(wf(ii));
        plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clr, 'LineStyle', 'none');
    end
end
