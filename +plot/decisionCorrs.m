function [pss, mdl] = decisionCorrs(vs, tp, ynm)
    figure; hold on;
    
    dts = unique({vs.dt});
    vd = vs(~[vs.isCell]);
    vs = vs(strcmp({vs.type}, tp));
    pss = [];
    for ii = 1:numel(dts)        
        vdc = vd(strcmp({vd.dt}, dts{ii}));
        w1 = vdc.w;
        vsc = vs(strcmp({vs.dt}, dts{ii}));
        xs = nan(1,numel(vsc));
        for jj = 1:numel(vsc)
            w2 = vsc(jj).w;
            if vsc(jj).targPref == 2
                w2 = -w2;
            end
            xs(jj) = corr(w1, w2);
        end
        ys = [vsc.(ynm)];
        [~,ix] = sort(xs);
        plot(xs(ix), ys(ix), 'color', [0.8 0.8 0.8]);
        pss = [pss; xs' ys'];
    end

    xs = pss(:,1); ys = pss(:,2);
    mdl = fitlm(xs, ys)
    plot(mdl, 'marker', 'o');
    legend off;
    
    vfcn = @(x) sprintf('%0.1f', x);
    nbins = 12;
    bins = linspace(floor(min(xs)),ceil(max(xs)),nbins);
    [~,ee] = histc(xs, bins);
    cents = bins(1:end-1) + diff(bins)/2;
    cents = cents(unique(ee));
%     boxplot(ys, ee, 'positions', cents, 'labels', ...
%         arrayfun(vfcn, cents, 'uni', 0));
    
    bs = linspace(floor(min(xs)),ceil(max(xs)), nbins-1);
    set(gca, 'xtick', bs);
    set(gca, 'xticklabel', arrayfun(vfcn, bs, 'uni', 0));
    xlim([floor(min(xs)) ceil(max(xs))]);
    xl = xlim; xlim(xl*1.1);
    ylim([floor(min(ys)) ceil(max(ys))]);

    xlabel(['corr(' tp ' RF, decision RF)']);
    ylabel(ynm);
    title([ynm ' vs. correlation between decision RF and ' tp]);
end
