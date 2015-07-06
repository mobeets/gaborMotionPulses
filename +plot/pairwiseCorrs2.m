function pairwiseCorrs2(vs, tp, xnm, ynm)
    vs = vs(strcmp({vs.type}, tp));
    dts = unique({vs.dt});
    ys = [];
    xs = [];
%     xnm = 'w';
%     ynm = 'Yfrz';
    for ii = 1:numel(dts)
        v0 = vs(strcmp({vs.dt}, dts{ii}));
        if numel(v0) < 2
            continue;
        end
        for jj = 1:numel(v0)
            for kk = 1:numel(v0)
                if jj == kk
                    continue;
                end
                
                Y1 = v0(jj).(ynm);
                Y2 = v0(kk).(ynm);
                
%                 dps = unique(v0(jj).dirprob);
%                 Y1h = nan(numel(dps),1);
%                 Y2h = nan(numel(dps),1);
%                 for ll = 1:numel(dps)
%                     Y1h(ll) = mean(Y1(v0(jj).dirprob == dps(ll)));
%                     Y2h(ll) = mean(Y2(v0(kk).dirprob == dps(ll)));
%                 end
%                 Y1 = Y1h;
%                 Y2 = Y2h;
                
                iy = ~isnan(Y1) & ~isnan(Y2);
                if sum(iy) == 0
                    continue;
                end

                X1 = v0(jj).(xnm);
                X2 = v0(kk).(xnm);
                ix = ~isnan(X1) & ~isnan(X2);
                if sum(ix) == 0
                    continue;
                end
                xs = [xs corr(X1(ix), X2(ix))];
                ys = [ys tools.dcor(Y1(iy), Y2(iy), false)];
            end
        end
    end
    figure;
    subplot(2,2,4); hold on;
    hist(xs, linspace(-1,1,21));
    subplot(2,2,1); hold on;
    hist(ys, linspace(0,1,10));
    subplot(2,2,2); hold on;
    scatter(abs(xs), ys);
    xlim([-1 1]);
    ylim([0 1]);

end
