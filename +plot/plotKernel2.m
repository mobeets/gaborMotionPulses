function plotKernel2(wf, xy, sz, fp, t1, t2, name)
    if nargin < 7
        name = '';
    end
    if nargin < 6 || any(isnan(t2))
        t2 = [nan nan];
    end
    if nargin < 5 || any(isnan(t1))
        t1 = [nan nan];
    end
    if nargin < 4 || any(isnan(fp))
        fp = [nan nan];
    end
    if nargin < 3 || isnan(sz)
        sz = 50;
    end
    [nw, nt] = size(wf);
    
    figure;
    clrs = plot.getColors(wf(:));
    hold on;
    axis off;
    set(gcf, 'color', 'w');    
    
    xy0 = xy;
    xshift = 0;
    shift = max(xy) - min(xy);
    for ii = 1:nt
        for jj = 1:nw
            clr = clrs((ii-1)*nw+jj,:);
            plot(xy(jj,1), xy(jj,2), 'Marker', '.', 'MarkerSize', sz, ...
                'Color', clr, 'LineStyle', 'none');
        end
        text(median(xy(:,1))-0.5, min(xy(:,2)) - 2, ['t_' num2str(ii)], ...
            'FontSize', 18);
        if ii < nt
            xy(:,1) = xy(:,1) + 1.5*shift(1);
        end
        if ii == 1
            if ~any(isnan(fp)) || ~any(isnan(t1)) || ~any(isnan(t2))
                mx = nanmax([t1(1) t2(1) fp(1)]);
                xshift = mx - max(xy0(:,1));
                if xshift > 0
%                     xshift = xshift - 2;
                    xy(:,1) = xy(:,1) + xshift;
                else
                    xshift = 0;
                end
            end
        end
    end
    
    if ~any(isnan(fp))
        scatter(fp(1), fp(2), sz, [0.8 0.2 0.2], 's', 'filled');
    end
    if ~any(isnan(t1))
        scatter(t1(1), t1(2), sz, [0.2 0.8 0.2], 'filled');
    end
    if ~any(isnan(t2))
        scatter(t2(1), t2(2), sz, [0.2 0.5 0.2], 'filled');    
    end

    title(name, 'FontSize', 18);
    set(gcf, 'Position', [100 100 1200 260]);
    axis equal;
end
