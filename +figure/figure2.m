%% Single Cell Prediction Comparison

exampleCellMT = '20140307-MT_3';
vMT = vuMT(strcmp({vuMT.name}, exampleCellMT));
d = io.loadDataByDate(vMT.dt, vMT.isNancy);
goodCells = find([vuMT.dPrime] > 1);

h = plot.trueVsModelSpikes(vMT);

figure(h(1))
idx = 300+[0 100];
xlim(idx)
set(h(1), 'PaperSize', [4 1], 'PaperPosition', [0 0 4 1])
title(strrep(vMT.name, '_', ' '))
saveas(h(1), fullfile(figDir, sprintf('trialCompare%d.pdf', exampleCellMT)))

figure(h(3))
set(h(3), 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
saveas(h(3), fullfile(figDir, sprintf('MotionDirSensitivity%d.pdf', exampleCellMT)))

%%
figure(110); clf
% tr = find(d.stim.dirprob==0,1);
tr = randi(d.stim.nTrials);
imagesc(squeeze(d.stim.pulses(tr,:,:))'); %colormap gray
xlabel('Pulse #')
ylabel('Gabor Id')
set(gcf, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
title(strrep(vMT.name, '_', ' '))
% saveas(gcf, fullfile(figDir, sprintf('stimTrial%d.pdf', kNeuron)))

%%
figure(111); clf
imagesc(reshape(d.stim.pulses, d.stim.nTrials, [])'); colormap gray
xlim(idx)
set(gcf, 'PaperSize', [4 1], 'PaperPosition', [0 0 4 1])
title(strrep(vMT.name, '_', ' '))
saveas(gcf, fullfile(figDir, sprintf('stimulus%d.pdf', exampleCellMT)))

%% r-sq before and after AR-2 model

figure;
set(gca, 'FontSize', 14);
set(gcf, 'color', 'w');
plot.fitSummaries(vuMT, 'score', 'score_AR');
