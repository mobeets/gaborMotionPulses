function mdl = boxScatterFitPlot(xs, ys, showBox, showFit, nbins, PCregression)
    if nargin < 5
        nbins = 12;
    end
    if nargin < 4
        showFit = true;
    end
    if nargin <6
        PCregression=true;
    end
    if nargin < 3
        showBox = true;
    end
    set(gca, 'FontSize', 14);
    hold on;
    
    % plot data
    scatter(xs, ys, 20, [0.8 0.8 0.8], 'filled');
    
    % smoothed running mean
%     plot(fit(xs, ys, 'poly4'), '-');
%     [~,idx] = sort(xs);
%     plot(xs(idx), smooth(xs(idx), ys(idx), 0.3, 'loess'), 'r-');
    
    % plot model fit
    mdl = fitlm(xs, ys)
    rsq = mdl.Rsquared.Ordinary;
    pval = mdl.Coefficients.pValue(2);
    ix = ~isnan(xs) & ~isnan(ys);
    rho = corr(xs(ix), ys(ix));
    cc = mdl.coefCI;
    cc = cc(2,:);
    if showFit
        hs = mdl.plot();
        set(hs(1), 'visible', 'off'); % hide data points
        set(hs(2:end), 'color', 'k');
        set(hs(2), 'LineWidth', 2);        
    end
    legend off;
    if PCregression
        X=[xs(:) ys(:)];
        xd=xlim;
        yd=ylim;
        coeff = pca(X);
        meanX=mean(X);
        Xfit = @(x,k) (coeff(1,k)/coeff(2,k))*(x-meanX(1))+meanX(2);
        cmap=lines;
        plot(xlim,Xfit(xlim,1), 'Color', cmap(1,:))
        plot(xlim,Xfit(xlim,2), 'Color', cmap(2,:))
        xlim(xd)
        ylim(yd)
    end
    
    title(['corr(xs,ys) = ' num2str(rho) ...
        ', fit r^2 = ' num2str(rsq) ', pval = ' num2str(pval) ...
        ', slope CI = ' num2str(cc)]);

    % plot box-whiskers
    vfcn = @(x) sprintf('%0.1f', x);
    minx = min(xs)-0.01;
    maxx = max(xs)+0.01;
    bins = linspace(minx, maxx, nbins);
    [~,ee] = histc(xs, bins);
    cents = bins(1:end-1) + diff(bins)/2;
    cents = cents(unique(ee));
    if showBox
        hs = boxplot(ys, ee, 'positions', cents, 'labels', ...
            arrayfun(vfcn, cents, 'uni', 0));
        set(hs(5,:), 'color', 'k');
        set(hs(6,:), 'color', 'r');
        set(hs(6,:), 'LineWidth', 2);
    end
    
    bs = linspace(floor(min(xs)),ceil(max(xs)), nbins-1);
    set(gca, 'xtick', bs);
    set(gca, 'xticklabel', arrayfun(vfcn, bs, 'uni', 0));
%     xlim([floor(min(xs)) ceil(max(xs))]);
    xl = xlim; xlim(xl*1.1);
%     ylim([floor(min(ys)) ceil(max(ys))]);
    
    box off;
    
end
