%% cartoon of color projection

plot.hyperFlowColorCartoon;
% plot.saveFigure('exampleOverlayProjection', figDir, gcf, 'png');

%% example hyperflow with arrows and no contourf

exampleCells = {'20150324a-MT_14', '20150324a-MT_13', ...
    '20140304-MT_2', '20140307-MT_1', '20140304-MT_3', '20150324a-MT_5'};

for kCell = 1:numel(exampleCells)
    exampleCellMT = exampleCells{kCell};
    vMT = vuMT(strcmp({vuMT.name}, exampleCellMT));
    if isempty(vMT)
        continue;
    end
    d = io.loadDataByDate(vMT.dt, vMT.isNancy);
    n = d.neurons{vMT.cellind};
    if isempty(n.hyperflow)
        continue
    end
    t1 = d.stim.targ1XY; t2 = d.stim.targ2XY;
    
    figure;
    ax = tight_subplot(1,1,.1,.1, .1);
    axes(ax) %#ok<LAXES>
    set(gcf, 'color', 'w');
    plot.getColors([0 1]);
    plot.plotHyperflowMT(n, t1, t2, false); axis xy
    xd = xlim;
    yd = ylim;
    plot(0,0, '+w', 'MarkerSize', 10, 'Linewidth', 2) % plot fixation cross
    xlim(xd) % rescale axes to the image size from before
    ylim(yd)
    title([n.exname '-' num2str(n.id, '%02.0f')]) % name of neuron (real id, not index)
    set(gca, 'XTick', [xd(1) 0 xd(end)], 'XTickLabel', round([xd(1) 0 xd(end)]), ...
        'YTick', [yd(1) 0 yd(end)], 'YTickLabel', round([yd(1) 0 yd(end)]))
    figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [50 50])    
%     plot.saveFigure(sprintf('exampleOverlay%s', exampleCellMT), figDir, gcf, 'png');
    
end

%% Plot Contour overlays and PSTHs
showTargs = false;
showHyperflow = true;
contourNoQuiver = true;

for ii = 1:numel(exampleCells)
    vMT = vuMT(strcmp({vuMT.name}, exampleCells{ii}));
    if numel(vMT) ~= 1
        continue;
    end
    
    % hyperflow ASD overlay
    d = io.loadDataByDate(vMT.dt, vMT.isNancy);
    n = d.neurons{vMT.cellind};
    figure;
    plot.plotSaccadeKernelOverlay(d.stim, n, vMT, showTargs, ...
        showHyperflow, contourNoQuiver);
    xd = xlim;
    yd = ylim;
    plot(0,0, '+k', 'MarkerSize', 10.5, 'Linewidth', 2.5) % plot fixation cross
    plot(0,0, '+w', 'MarkerSize', 10, 'Linewidth', 2) % plot fixation cross
    xlim(xd) % rescale axes to the image size from before
    ylim(yd)
    title([n.exname '-' num2str(n.id, '%02.0f')]) % name of neuron (real id, not index)
    set(gca, 'XTick', [xd(1) 0 xd(end)], 'XTickLabel', round([xd(1) 0 xd(end)]), ...
        'YTick', [yd(1) 0 yd(end)], 'YTickLabel', round([yd(1) 0 yd(end)]))
    figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [50 50])
    plot.saveFigure(sprintf('hyperflowRF%s', vMT.name), figDir, gcf, 'png');


    % PSTH split by motion direction
%     figure; hold on;
%     set(gcf, 'color', 'w');
%     set(gca, 'FontSize', 14);
%     targPref = nanmax([n.targPref, 1]);
%     event = sum(sum(d.stim.pulses, 3), 2) > 0; event = -(event-2);
%     plot.psthByEvent(d, n, event == targPref, ii, {'-mot', '+mot'});
%     ylabel('spikes/sec by motion dir');
            
%     plot.saveFigure(['psthByMot_' vMT.name], figDir, gcf, 'pdf');
end
