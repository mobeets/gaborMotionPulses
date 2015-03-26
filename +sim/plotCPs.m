function plotCPs(vals, xs, xslbl)
    marker1 = '-';
    marker2 = '--';
    [~, idx] = sort(xs);
    if isempty(xslbl)
        xs = 1:max(idx);
    else
        xs = xs(idx);
    end
    
    plotVal(xs, vals, idx, 'cp_Ypos', marker1, 2, 'b', 'Y_{pos}');
    plotVal(xs, vals, idx, 'cp_Yneg', marker1, 2, 'r', 'Y_{neg}');
    plotVal(xs, vals, idx, 'cp_Yh', marker1, 2, [0.3 0.3 0.3], 'Y_{pos} + Y_{neg}');
    plotVal(xs, vals, idx, 'cp_Y', marker1, 3, 'k', 'Y_{pos} + Y_{neg} + Y_{res}');
    plotVal(xs, vals, idx, 'cp_Yposres', marker2, 2, [0.4 0.4 1], 'Y_{pos} + Y_{res}');
    plotVal(xs, vals, idx, 'cp_Ynegres', marker2, 2, [1 0.4 0.4], 'Y_{neg} + Y_{res}');
    plotVal(xs, vals, idx, 'cp_Yres', marker2, 2, [0.7 0.7 0.7], 'Y_{res}');
    
%     plotVal(xs, vals, idx, 'cp_Yadj1', marker1, 2, 'g', 'Y_{pos} - Y_{neg}');
%     plotVal(xs, vals, idx, 'cp_Yadj2', marker2, 2, 'g', 'Y_{pos} - Y_{neg} + Y_{res}');
    if ~isempty(xslbl)
        xlabel(xslbl);
    else
        xlabel('index');
    end
    ylabel('CP');
    ylim([0 1]);
    legend('Location', 'NorthEastOutside');
end

function plotVal(xs, vals, idx, name, marker, lw, clr, lbl)
    ys = [vals(idx).(name)];
    plot(xs, ys, marker, 'LineWidth', lw, 'Color', clr, ...
        'DisplayName', ['CP(' lbl ')']);
end
