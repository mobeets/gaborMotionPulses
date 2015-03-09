function fig = plotKernel(xy, wf, vmax, figLbl, sz, figSz, clrFcn, xLblFcn, yLblFcn)
% plots an nw-by-nt spatiotemporal kernel
%   creates nt subplots each with nw weights
% 
% xy - spatial coords or wf
% wf - weights to plot
% vmax - normalizer for wf (default is maximum value in wf)
% sz - size of markers
% figSz - size of figure
% figLbl - prefix for title of each subplot
% clrFcn - color function handle of form color = f(wf(i,j)) for any i,j
% lblFcn - generates labels for each subplot
% 
    assert(~any(isnan(wf(:))));
    if nargin < 9 || ~isa(yLblFcn, 'function_handle')
        yLblFcn = @(ii) '';
    end
    if nargin < 8 || ~isa(xLblFcn, 'function_handle')
        xLblFcn = @(ii) ['t=' num2str(ii)];
    end
    if nargin < 7 || ~isa(clrFcn, 'function_handle')
        clrFcn = plot.colorScheme();
    end
    if nargin < 6 || isnan(figSz)
        figSz = 1.0;
    end
    if nargin < 5 || isnan(sz)
        sz = 50;
    end
    if nargin < 4 || any(isnan(figLbl))
        figLbl = '';
    end
    if nargin < 3 || isnan(vmax)
        vmax = max(abs(wf(:)));
    end    
    mrg = 1.0;
    titleFontSize = 14;
    
    [nw, nt] = size(wf);
    fig = figure;
    ha = plot.tight_subplot(1, nt, [.01 .03], [.1 .01], [.01 .01]);

    wf = wf/vmax; % normalize
    for ii = 1:nt
        axes(ha(ii)); hold on;
        for jj = 1:nw
            clr = clrFcn(wf(jj,ii));
            plot(xy(jj,1), xy(jj,2), 'Marker', '.', 'MarkerSize', sz, ...
                'Color', clr, 'LineStyle', 'none');
        end
        subplotFormat();
        xlim([min(xy(:,1))-mrg, max(xy(:,1))+mrg]);
        ylim([min(xy(:,2))-mrg, max(xy(:,2))+mrg]);
        if ii == round(nt/2)
            ht = title(figLbl);
            set(ht, 'FontSize', titleFontSize);
        end
        xlabel(xLblFcn(ii)); % acts like a subplot title
        ylabel(yLblFcn(ii));
    end
    plot.suplabel(figLbl, 't'); 
    
    pos = get(gcf,'Position');
    pos(3:4) = figSz*[10e2 2e2];
    set(gcf, 'Position', pos);
    set(gcf,'color','w');
end

function subplotFormat()
    axis square;
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'LineWidth', 1)
    box on;
end

