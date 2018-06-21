function fig = showAllRfs(cells)

    n = numel(cells);
    ncols = ceil(sqrt(n));
    nrows = ceil(n/ncols);
    nrows = 12; ncols = 11;
    offset = 5;    
    sz = 50;
    cmin = nan;
    
    smoothness = 1./[cells.rfSpatialVariability];
    smoothness = sign([cells.dPrime]).*smoothness;
    [~,ix] = sort(smoothness);
    cells = cells(ix);
    
    fig = plot.init;
    c = 0;
    for jj = 1:nrows
        for ii = 1:ncols
            c = c + 1;
            if c > numel(cells)
                break;
            end
            
            cell = cells(c);
            rf = cell.wsep.spatial_RF;
            if sum(cell.wsep.scalar_RF*cell.wsep.temporal_RF) < 0
                rf = -rf;
            end
            clrs = plot.getColors(rf, true, cmin);
            
            xy = cell.Xxy;
            xy = zscore(xy);
            cx = (ii-1)*offset;
            cy = -(jj-1)*offset;
            scatter(xy(:,1) + cx, xy(:,2) + cy, sz, clrs, 'filled');
        end
    end    
    axis off;
    axis equal;
    plot.setPrintSize(gcf, struct('width', 6, 'height', 6));

end
