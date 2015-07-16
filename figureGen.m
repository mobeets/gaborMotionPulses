%% load fit data

fitdir = '20150615';
fitbasedir = 'data';
figDir = '~/Dropbox/gaborMotionPulseASD/figures/Figure02_ASDmethods/';
ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
vu = [vp vn];

%% pairwise plots

Yh = 'Ylow';
% vut = vu([vu.nzertrials] >= 30);

% noise-correlation vs. distance of RF centers
plot.pairwiseCorrs(vu([vu.isMT]), 'rf_center', Yh, ...
    nan,nan,nan,@(x,y) norm(x-y,2));
% plot.pairwiseCorrs(vut([vut.isMT]), 'rf_center', Yh, ...
%     nan,nan,nan,@(x,y) norm(x-y,2));
set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
plot.saveFigure('MT - noise-corr vs. rf-center', 'tmp', gcf, 'pdf');



% noise-correlation vs. RF similarity
% plot.pairwiseCorrs(vu([vu.isLIP]), 'w', Yh);
% plot.saveFigure('LIP - noise-corr vs. rf-corr', 'tmp', gcf);
plot.pairwiseCorrs(vu([vu.isMT]), 'w', Yh);
% plot.pairwiseCorrs(vut([vut.isMT]), 'w', Yh);
set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
plot.saveFigure('MT - noise-corr vs. rf-corr', 'tmp', gcf, 'pdf');

%% CP plots

cpY = 'cp_Ylow';
% vut = vu([vu.nzertrials] >= 30);
vut = vu(arrayfun(@(v) all(v.nlowmottrials > 30), vu));
% vut = vut(round([vut.(cpY)]) ~= [vut.(cpY)]);

% plot.boxScatterFitPlotWrap(vu([vu.isMT]), 'dPrime', cpY);
plot.boxScatterFitPlotWrap(vut([vut.isMT]), 'dPrime', cpY);
plot.saveFigure(['MT - ' cpY ' vs. dPrime'], 'tmp', gcf);

plot.boxScatterFitPlotWrap(vut([vut.isMT]), 'rf_ecc', cpY, ...
    true, true, 6);
% xlim([0 0.6]);
plot.saveFigure(['MT - ' cpY ' vs. rf-eccentricity'], 'tmp', gcf);

% plot.boxScatterFitPlotWrap(vu([vu.isMT]), 'dPrime', 'rf_ecc', ...
%     false, false, 6);

%%

% far away
% low noise corr
% high d'

vut = vu([vu.nzertrials] >= 30);
plot.pairwiseCorrs(vut([vut.isMT]), 'rf_center', 'Yfrz', ...
    nan,nan,nan,@(x,y) norm(x-y,2));

%% Single Cell Prediction Comparison
goodCells = find([vu.dPrime] > 1 & [vu.isMT]);
k = 6;
kNeuron = goodCells(k);
close all
h = figure.trueVsModelSpikes(vu(kNeuron));

figure(h(1))
idx = 300+[0 100];
xlim(idx)
set(h(1), 'PaperSize', [4 1], 'PaperPosition', [0 0 4 1])
title(strrep(vu(kNeuron).name, '_', ' '))
saveas(h(1), fullfile(figDir, sprintf('trialCompare%d.pdf', kNeuron)))

figure(h(3))
set(h(3), 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
saveas(h(3), fullfile(figDir, sprintf('MotionDirSensitivity%d.pdf', kNeuron)))

%% get stimulus data
d = io.loadDataByDate(vu(kNeuron).dt, vu(kNeuron).isNancy);
%%
figure(110); clf
% tr = find(d.stim.dirprob==0,1);
tr = randi(d.stim.nTrials);
imagesc(squeeze(d.stim.pulses(tr,:,:))'); %colormap gray
xlabel('Pulse #')
ylabel('Gabor Id')
set(gcf, 'PaperSize', [4 4], 'PaperPosition', [0 0 4 4])
title(strrep(vu(kNeuron).name, '_', ' '))
% saveas(gcf, fullfile(figDir, sprintf('stimTrial%d.pdf', kNeuron)))
%%
figure(111); clf
imagesc(reshape(d.stim.pulses, d.stim.nTrials, [])'); colormap gray
xlim(idx)
set(gcf, 'PaperSize', [4 1], 'PaperPosition', [0 0 4 1])
title(strrep(vu(kNeuron).name, '_', ' '))
saveas(gcf, fullfile(figDir, sprintf('stimulus%d.pdf', kNeuron)))