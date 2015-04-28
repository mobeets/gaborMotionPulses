function noiseCorrByStim(v1, v2, xs, xstr, ystr, nxbins)
% example: plot.histSpikesByStim(v.dirprob, v.Y)
% 

    X1 = v1.(xstr);
    X2 = v2.(xstr);
    
    if nargin > 5
        xs = linspace(min([X1; X2]), max([X1; X2]), nxbins);
        [~, X1] = histc(X1, xs);
        [~, X2] = histc(X2, xs);
        X1 = xs(X1);
        X2 = xs(X2);
    end
    
    Y1 = v1.(ystr);
    Y2 = v2.(ystr);

    ys = nan(numel(xs), 1);
    for ii = 1:numel(xs)
        ys(ii) = corr(Y1(X1 == xs(ii)), Y2(X2 == xs(ii)));
    end
    figure; hold on;
    plot([min(xs) max(xs)], [0 0], 'k--');
    scatter(xs, ys, 'filled');
    xlabel('X');
    ylabel('noise correlation');
end
