function [pss, mdl] = corrs(vs, xnm, ynm, lbl, xlbl, ylbl)
    if nargin < 6
        ylbl = ynm;
    end
    if nargin < 5
        xlbl = xnm;
    end
    if nargin < 4
        lbl = '';
    end
    lblfcn = @(x) ['corr(' x '_1, ' x '_2)'];
    xlbl = lblfcn(xlbl);
    ylbl = lblfcn(ylbl);

    dts = unique({vs.dt});
    pss = [];
    figure;
    subplot(3,1,1); hold on;
    set(gca, 'FontSize', 14);
    for ii = 1:numel(dts)
        vts = vs(strcmp({vs.dt}, dts{ii}));
        ps = [];
        for jj = 1:numel(vts)
            if jj+1 > numel(vts)
                continue;
            end
            v0 = vts(jj);
            for kk = jj+1:numel(vts)
                v1 = vts(kk);

                xs = v0.(xnm);
                ys = v1.(xnm);
%                 xs = v0.wfSvd_2(:);
%                 ys = v1.wfSvd_2(:);
%                 xs = v0.wfSvd_U(:,1);
%                 ys = v1.wfSvd_U(:,1);

                ix = ~isnan(xs) & ~isnan(ys);
                xs = xs(ix);
                ys = ys(ix);
                if numel(xs) == 0
                    continue;
                end
                p1 = corr(xs, ys);

                xs = v0.(ynm);
                ys = v1.(ynm);
                ix = ~isnan(xs) & ~isnan(ys);
                xs = xs(ix);
                ys = ys(ix);
                if numel(xs) == 0
                    continue;
                end
                p2 = corr(xs, ys);
                ps = [ps; p1 p2];
            end    
        end
        if numel(ps) > 1
            [~,ix0] = sort(ps(:,1));
            plot(ps(ix0,1), ps(ix0,2), 'Color', [0.8 0.8 0.8]);
            pss = [pss; ps];
        end
    end

    xs = pss(:,1); ys = pss(:,2);
    mdl = fitlm(xs, ys)
    plot(mdl, 'marker', '.');
    legend off;
    
    vfcn = @(x) sprintf('%0.1f', x);
    nbins = 12;
    bins = linspace(floor(min(xs)),ceil(max(xs)),nbins);
    [~,ee] = histc(xs, bins);
    cents = bins(1:end-1) + diff(bins)/2;
    cents = cents(unique(ee));
    boxplot(ys, ee, 'positions', cents, 'labels', ...
        arrayfun(vfcn, cents, 'uni', 0));
    
    bs = linspace(floor(min(xs)),ceil(max(xs)), nbins-1);
    set(gca, 'xtick', bs);
    set(gca, 'xticklabel', arrayfun(vfcn, bs, 'uni', 0));
    xlim([floor(min(xs)) ceil(max(xs))]);
    xl = xlim; xlim(xl*1.1);
    ylim([floor(min(ys)) ceil(max(ys))]);

    xlabel(xlbl);
    ylabel(ylbl);
    
%     vfcn = @(x) sprintf('%0.2f', x);
%     lbl = ['mean(' xlbl ') = ' vfcn(nanmean(xs)) ...
%         ', mean(' ylbl ') = ' vfcn(nanmean(ys))];
    title(lbl);
    
    subplot(3,1,2); hold on;
    bins = linspace(floor(min(xs)),ceil(max(xs)),nbins+1);
    hist(xs, bins);
    xlabel(xlbl);
    xlim(1.1*[min(bins) max(bins)]);
    subplot(3,1,3); hold on;
    bins = linspace(floor(min(ys)),ceil(max(ys)),nbins+1);
    hist(ys, bins);
    xlim(1.1*[min(bins) max(bins)]);
    xlabel(ylbl);
    
%     subplot(3,1,1); hold on;
end
