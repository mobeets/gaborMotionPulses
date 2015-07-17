%% cartoon of color projection

plot.hyperFlowColorCartoon;
plot.saveFigure('exampleOverlayProjection', figDir, gcf, 'pdf');

%% example hyperflow with arrows and no contourf

% exampleCellMT = '20140304-MT_3';
exampleCells = [cellstr(num2str((1:7)', '20140304-MT_%d')); cellstr(num2str((1:2)', '20140307-MT_%d')); cellstr(num2str((1:10)', '20150324a-MT_%d'))];
for kCell = 1:numel(exampleCells)
    exampleCellMT = exampleCells{kCell};
    vMT = vuMT(strcmp({vuMT.name}, exampleCellMT));
    if isempty(vMT)
        continue
    end
    d = io.loadDataByDate(vMT.dt, vMT.isNancy);
    n = d.neurons{vMT.cellind};
    if isempty(n.hyperflow)
        continue
    end
    t1 = d.stim.targ1XY; t2 = d.stim.targ2XY;
    figure;
    hold on;
    set(gcf, 'color', 'w');
    axis off;
    plot.getColors([0 1]);
    plot.plotHyperflowMT(n, t1, t2, false);
%     scatter(t1(:,1), t1(:,2), 50, [0.2 0.8 0.2], 'filled');
%     scatter(t2(:,1), t2(:,2), 50, [0.2 0.5 0.2], 'filled');
%     axis equal;
    
    plot.saveFigure(sprintf('exampleOverlay%s', exampleCellMT), figDir, gcf, 'pdf');
end

%%

showTargs = false;
showHyperflow = true;
contourNoQuiver = true;

% exampleCells = cellstr(num2str((1:7)', '20140304-MT_%d'));
exampleCells = [cellstr(num2str((1:7)', '20140304-MT_%d')); cellstr(num2str((1:2)', '20140307-MT_%d')); cellstr(num2str((1:30)', '20150324a-MT_%d'))];
% exampleCells = {'20140305-MT_1', '20140218-MT_1', '20140304-MT_3', ...
%     '20150324a-MT_4', '20150324a-MT_5', '20150324a-MT_6'};
% exampleCells = {'20150519-MT_5', '20150519-MT_6'};
% exampleCells = {'20150304b-MT_14', '20150304b-MT_8'};

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
    axis off;
    plot.saveFigure(['hyperflowRF_' vMT.name], figDir, gcf, 'pdf');

    % PSTH split by motion direction
    figure; hold on;
    set(gcf, 'color', 'w');
    set(gca, 'FontSize', 14);
    targPref = nanmax([n.targPref, 1]);
    event = sum(sum(d.stim.pulses, 3), 2) > 0; event = -(event-2);
    plot.psthByEvent(d, n, event == targPref, ii, {'-mot', '+mot'});
    ylabel('spikes/sec by motion dir');
            
    plot.saveFigure(['psthByMot_' vMT.name], figDir, gcf, 'pdf');
end
