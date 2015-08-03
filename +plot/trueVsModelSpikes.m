function [h,S] = trueVsModelSpikes(v, idx)
    if ~exist('idx', 'var')
        idx = 1:100;
    end
    Y = v.Y;
    Yh = v.Yh;
    Yh1 = v.YhAR;
    ix = ~isnan(Y); Y = Y(ix); Yh = Yh(ix); Yh1 = Yh1(ix);
    
    Fs = 5;
    h(1) = figure;
    ax = tight_subplot(1,1,.1, .1, .1);
    axes(ax)
    % output for pyplot
    S = struct();
    Yp = Y(idx);
    Yhp = Yh1(idx);
    xx = [1 1:numel(Yp) numel(Yp)];
    yy = [0 Yp' 0];
    fill(xx(:),yy(:),'k', 'FaceColor', .5*[1 1 1], 'Linestyle', 'none'); hold on
    plot(Yhp, 'r', 'Linewidth', 1);   
    set(gcf, 'color', 'w');
    xlabel('trial'); ylabel('spike count'); title(v.name);
    
    S.trial_prediction.Y = Y;
    S.trial_prediction.Yhat = Yh1;
    S.trial_prediction.ix   = idx;
    S.trial_prediction.fillx = xx;
    S.trial_prediction.filly = yy;
    S.trial_prediction.xstim = v.dirprob(ix);
        
    xlim([0 numel(Yp)])
    h(2) = figure;
    ax = tight_subplot(1,1,.1, .1, .1);
	set(gcf, 'color', 'w'); 
    axes(ax)
    set(gca, 'FontSize', Fs);
    xlabel('dirprob'); ylabel('spike count'); title(v.name);
    xs = v.dirprob(ix);
    scatter(xs, Y, 'b');hold on;
    scatter(xs, Yh1, 'r');
    
    
    h(3) = figure; 
    ax = tight_subplot(1,1,1, .1, .1);
    hold on; set(gcf, 'color', 'w');
    axes(ax)
    set(gca, 'FontSize', Fs);
    xlabel('dirstrength'); ylabel('spike count'); title(v.name);
    xs = v.dirstrength(ix);  
    [xso, ix] = sort(xs);
    plot.boxScatterFitPlot(xs, Y, true, false, 10)
    plot(xso, smooth(xso, Yh1(ix), 0.3, 'loess'), 'b-');
end
