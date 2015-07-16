function h = trueVsModelSpikes(v)
    Y = v.Y;
    Yh = v.Yh;
    Yh1 = v.YhAR;
    ix = ~isnan(Y); Y = Y(ix); Yh = Yh(ix); Yh1 = Yh1(ix);
    
    Fs = 5;
    h(1) = figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', Fs);
    xlabel('trial'); ylabel('spike count'); title(v.name);
    xx = [1 1:numel(Y) numel(Y)];
    yy = [0 Y' 0];
    plot(xx,yy,'k');
%     fill(xx,yy,'k', 'FaceColor', .5*[1 1 1])
%     stairs(Y, 'k');
%     stairs(Yh, 'b');
    plot(Yh1, 'r', 'Linewidth', 1);    

    h(2) = figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', Fs);
    xlabel('dirprob'); ylabel('spike count'); title(v.name);
    xs = v.dirprob(ix);
    scatter(xs, Y, 'b');
    scatter(xs, Yh1, 'r');
    

    h(3) = figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', Fs);
    xlabel('dirstrength'); ylabel('spike count'); title(v.name);
    xs = v.dirstrength(ix);  
    [xso, ix] = sort(xs);
    plot.boxScatterFitPlot(xs, Y, true, false, 10)
    plot(xso, smooth(xso, Yh1(ix), 0.3, 'loess'), 'b-');
end
