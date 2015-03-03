function summaryByCell(dt, cellind, fitdir)
    [vals, data] = io.summariesByDate(dt, fitdir, 1);
    neuron = data.neurons{cellind};
    event = data.stim.targchosen;
    [Z, bins] = io.psthByEvent(data.stim, neuron, event);
    
    ind = strcmp({vals.name}, ['cell_' num2str(cellind)]);
    wf = vals(ind).mu;
    [u,s,v] = svd(wf);
    mu = u(:,1)*s(1)*v(:,1)';
    sc = vals(ind).score;
    sep = vals(ind).separability;
    fmt = @(val) sprintf('%0.2f', val);
    
    figure;
    set(gcf,'color','w');
    
    subplot(2, 2, 1); hold on;
    lw = 1;
    for ii = 1:numel(Z)
        if ii == 1
            lbl = 'pref';
        else
            lbl = 'anti';
        end
        plot(bins(:,1)*1000, Z{ii}, repmat('-', 1, ii), ...
            'Color', 'k', 'LineWidth', ii*lw, ...
            'DisplayName', lbl);
    end
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
    
    subplot(2, 2, 3); hold on;
    vmax = max(abs(wf(:)));
    plotKernelMini(data.Xxy, mu(:,1), vmax);
    title('spatial weights');
    xlabel(['sc=' fmt(sc) ', sep=' fmt(sep)]);
    
    subplot(2, 2, 4); hold on;
    bar(v(:,1), 'FaceColor', [0.3 0.8 0.3]);
    title('temporal weights');
    xlabel('pulse');

end

function plotKernelMini(xy, wf, vmax)
    sz = 50;
    clrFcn = plot.colorScheme();
    wf = wf/vmax;
    hold on;
    for ii = 1:numel(wf)
        clr = clrFcn(wf(ii));
        plot(xy(ii,1), xy(ii,2), 'Marker', '.', 'MarkerSize', sz, ...
            'Color', clr, 'LineStyle', 'none');
    end
end
