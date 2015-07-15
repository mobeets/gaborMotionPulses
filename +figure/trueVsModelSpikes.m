function trueVsModelSpikes(v)
    Y = v.Y;
    Yh = v.Yh;
    Yh1 = v.YhAR;
    ix = ~isnan(Y); Y = Y(ix); Yh = Yh(ix); Yh1 = Yh1(ix);
    
    figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', 14);
    xlabel('trial'); ylabel('spike count'); title(v.name);
    plot(Y, 'k');
    plot(Yh, 'b');
    plot(Yh1, 'r');    

    figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', 14);
    xlabel('dirprob'); ylabel('spike count'); title(v.name);
    xs = v.dirprob(ix);
    scatter(xs, Y, 'b');
    scatter(xs, Yh1, 'r');

    figure; hold on; set(gcf, 'color', 'w'); set(gca, 'FontSize', 14);
    xlabel('dirstrength'); ylabel('spike count'); title(v.name);
    xs = v.dirstrength(ix);
    scatter(xs, Y, 'b');
    scatter(xs, Yh1, 'r');
end
