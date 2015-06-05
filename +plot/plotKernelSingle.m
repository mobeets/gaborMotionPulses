function plotKernelSingle(xy, wf, vmax, sz)
    if nargin < 4 || isnan(sz)
        sz = 200;
    end
    if nargin < 3 || isnan(vmax)
        vmax = max(abs(wf(:)));
    end
    clrFcn = plot.colorScheme();
    wf = wf/vmax;
    hold on;
%     axis off;
    axis equal;
    set(gcf,'color','w');
    for ii = 1:numel(wf)
        clr = clrFcn(wf(ii));
        plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clr, 'LineStyle', 'none');
    end
end
