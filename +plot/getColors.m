function clrs = getColors(C, loadCmap, cmin, cmax)
    if nargin < 2
        loadCmap = true;        
    end
    if nargin < 3
        % diverging colorscheme, so take abs        
        cmax = max(abs(C));
        cmin = -cmax;
    end
    if loadCmap
        colormap(cbrewer('div', 'RdBu', 201, 'pchip'));
    end
    cmap = colormap;
    m = size(cmap,1);
    cind = fix((C-cmin)/(cmax-cmin)*m)+1;
    cind(cind<1) = 1;
    cind(cind>m) = m;
    clrs = cmap(cind,:);
%     caxis([cmin cmax]);
end
