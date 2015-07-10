function [pss, mdl] = pairwiseCorrs(vs, xnm, ynm, lbl, xlbl, ylbl, ...
    xfcn, yfcn)
% xfcn (yfcn) - how to compare two neurons for x-axis (y-axis)
%   - default is @corr
% 
    if nargin < 8
        yfcn = @corr;
    end
    if nargin < 7
        xfcn = @corr;
    end
    if nargin < 6 || all(isnan(ylbl))
        ylbl = ynm;
    end
    if nargin < 5 || all(isnan(xlbl))
        xlbl = xnm;
    end
    if nargin < 4 || all(isnan(lbl))
        lbl = '';
    end
    lblfcn = @(x) ['corr(' x '_1, ' x '_2)'];
    xlbl = lblfcn(xlbl);
    ylbl = lblfcn(ylbl);
    
    figure;
    set(gcf, 'color', 'w');
    subplot(3,1,1); hold on;

    dts = unique({vs.dt});
    pss = [];
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
%                 xs = v0.wfSvd_1(:);
%                 ys = v1.wfSvd_1(:);
%                 xs = v0.wfSvd_U(:,1);
%                 ys = v1.wfSvd_U(:,1);

                ix = ~isnan(xs) & ~isnan(ys);
                xs = xs(ix);
                ys = ys(ix);
                if numel(xs) == 0
                    continue;
                end
                p1 = xfcn(xs, ys);

                xs = v0.(ynm);
                ys = v1.(ynm);
                ix = ~isnan(xs) & ~isnan(ys);
                xs = xs(ix);
                ys = ys(ix);
                if numel(xs) == 0
                    continue;
                end
                p2 = yfcn(xs, ys);
                ps = [ps; p1 p2];
            end    
        end
        if numel(ps) > 1
            [~,ix0] = sort(ps(:,1));
            plot(ps(ix0,1), ps(ix0,2), 'Color', [0.8 0.8 0.8]);
            pss = [pss; ps];
        end
    end
    xs = pss(:,1);
    ys = pss(:,2);

    set(gca, 'FontSize', 14);
    mdl = plot.boxScatterFitPlot(pss(:,1), pss(:,2));
    xlabel(xlbl);
    ylabel(ylbl);
    
    subplot(3,1,2); hold on;
    nbins = 12;
    set(gca, 'FontSize', 14);
    bins = linspace(floor(min(xs)),ceil(max(xs)),nbins+1);
    hist(xs, bins);
    xlabel(xlbl);
    xlim(1.1*[min(bins) max(bins)]);
    
    subplot(3,1,3); hold on;
    set(gca, 'FontSize', 14);
    bins = linspace(floor(min(ys)),ceil(max(ys)),nbins+1);
    hist(ys, bins);
    xlim(1.1*[min(bins) max(bins)]);
    xlabel(ylbl);
    
end
