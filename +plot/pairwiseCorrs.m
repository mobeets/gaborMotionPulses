function [pss, mdl, nms] = pairwiseCorrs(vs, xnm, ynm, lbl, xlbl, ylbl, ...
    xfcn, yfcn, pairPredFcn)
% xfcn (yfcn) - how to compare two neurons for x-axis (y-axis)
%   - default is @corr
% pairPredFcn - given a pair of neurons, returns bool of whether to include
%   them or ignore
% 
    if nargin < 9
        pairPredFcn = @(v0, v1) true;
    end
    if nargin < 8
        yfcn = @corr;
    end
    if nargin < 7
        xfcn = @corr;
    end
    lblfcn = @(x) ['corr(' x '_1, ' x '_2)'];
    if nargin < 6 || all(isnan(ylbl))
        ylbl = lblfcn(ynm);
    end
    if nargin < 5 || all(isnan(xlbl))
        xlbl = lblfcn(xnm);
    end
    if nargin < 4 || all(isnan(lbl))
        lbl = '';
    end
    
    figure;
    set(gcf, 'color', 'w');
    subplot(3,1,1); hold on;

    dts = unique({vs.dt});
    pss = [];
    dps = [];
    nms = {};
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
%                 xs = v0.wfSvd_1(:,1);
%                 ys = v1.wfSvd_1(:,1);
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
%                 ps = [ps; p1 p2];
                
                if pairPredFcn(v0, v1)
                    ps = [ps; p1 p2];
%                     dps = [dps; p1 p2 v0.dPrime v1.dPrime v0.cp_Yfrz v1.cp_Yfrz];
                    nms = [nms; {v0.name, v1.name}];
                end
            end    
        end
        if numel(ps) > 1
            [~,ix0] = sort(ps(:,1));
%             plot(ps(ix0,1), ps(ix0,2), 'Color', [0.8 0.8 0.8]);
            pss = [pss; ps];
        end
    end
    
%     ix = dps(:,1) > nanmedian(dps(:,1)) & ... % rf-dist
%         dps(:,2) < nanmedian(dps(:,2)) & ... % noise-corr
%         dps(:,3) > nanmedian([vs.dPrime]) & ... % dprime-1
%         dps(:,4) > nanmedian([vs.dPrime]); % dprime-2
%     dps2 = dps(ix,:);
%     nms2 = nms(ix,:);
%     
    
    xs = pss(:,1);
    ys = pss(:,2);

    set(gca, 'FontSize', 14);
    mdl = plot.boxScatterFitPlot(pss(:,1), pss(:,2), false, false);
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
