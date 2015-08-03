%% Single Cell Prediction Comparison

exampleCellMT = '20140304-MT_12';
vMT = vuMT(strcmp({vuMT.name}, exampleCellMT));
d = io.loadDataByDate(vMT.dt, vMT.isNancy);
goodCells = find([vuMT.dPrime] > 1);

idx = 300:400;
[h, pyS] = plot.trueVsModelSpikes(vMT, idx);

% sample trials model fit quality
figure(h(1))
annotation('textbox', [.1 .8 .2 .05], 'String', strrep(vMT.name, '_', ' '), 'Linestyle', 'None')
title('')
figure.cleanupForPrint(h(1), 'FontSize', 8, 'PaperSize', [80 25])
saveas(h(1), fullfile(figDir, sprintf('trialCompare%d.png', exampleCellMT)))

% for saving out vector graphics with patch objects
% figure.cleanupForPrint(h(1), 'FontSize', 6, 'PaperSize', [80 25], true)
% plot2svg('~/Desktop/test2.svg', h(1))

% Spike Rate as a function of Net Motion
figure(h(3))
figure.cleanupForPrint(h(3), 'FontSize', 8, 'PaperSize', [25 25])
saveas(h(3), fullfile(figDir, sprintf('MotionDirSensitivity%s.pdf', exampleCellMT)))

%% Space-time image of a single trial
figure(110); clf
% tr = randi(d.stim.nTrials)
tr = 725;
imagesc(squeeze(d.stim.pulses(tr,:,:))', [-1 1]); colormap gray
set(gca, 'Xtick', [1 d.stim.nPulses], 'Ytick', [1 d.stim.nGabors])
xlabel('Pulse #')
ylabel('Gabor Id')
figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [25 25])
plot.saveFigure(sprintf('stimTrial%d.png', tr), figDir, gcf);

pyS.sample_trial.xyimage = squeeze(d.stim.pulses(tr,:,:));
pyS.sample_trial.gaborXY = vMT.Xxy;
%% All the trials
figure(111); clf
imagesc(reshape(d.stim.pulses(idx,:,:), numel(idx), [])'); colormap gray
set(gca, 'XTick', [1 numel(idx)], 'Ytick', [1 d.stim.nGabors*d.stim.nPulses])
xlabel('Trial')
ylabel('S(Trial)')
figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [80 25])
plot.saveFigure(sprintf('stimulus%d.png', exampleCellMT), figDir, gcf)

pyS.all_trials.image = reshape(d.stim.pulses, size(d.stim.pulses,1), [])';
%% r-sq before and after AR-2 model
figure(112); clf
set(gcf, 'color', 'w');
pyS.fit_summary= plot.fitSummaries(vuMT, 'b', 'score_AR', nan, false);
xlim([0 100])
hold on
xlabel('Mean Spike Rate')
ylabel('r^2')
figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [25 25])
plot.saveFigure('modelPerformancePopulation', figDir, gcf)

save(fullfile(figDir, 'figure02_pystruct.mat'), '-v7.3', '-struct', 'pyS')
