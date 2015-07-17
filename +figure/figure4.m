
Yh = 'YresAR';
cpY = 'cp_YresAR';

% vut = vuMT([vuMT.nzertrials] >= 30);
% vut = vuMT(arrayfun(@(v) all(v.nlowmottrials > 30), vuMT));
vut = vuMT;

%% noise-correlation vs. distance of RF centers
plot.pairwiseCorrs(vut, 'rf_center', Yh, ...
    nan,nan,nan,@(x,y) norm(x-y,2));
set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
plot.saveFigure('MT - noise-corr vs. rf-center', figDir, gcf, 'pdf');

%% noise-correlation vs. RF similarity
plot.pairwiseCorrs(vut, 'w', Yh);
set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
plot.saveFigure('MT - noise-corr vs. rf-corr', figDir, gcf, 'pdf');

%% CP vs. dPrime

plot.boxScatterFitPlotWrap(vut, 'dPrime', cpY);
plot.saveFigure(['MT - ' cpY ' vs. dPrime'], figDir, gcf);

%% CP vs. rfEcc
plot.boxScatterFitPlotWrap(vut, 'rf_ecc', cpY, ...
    true, true, 6);
plot.saveFigure(['MT - ' cpY ' vs. rf-eccentricity'], figDir, gcf);

%% dPrime vs. rfEcc
plot.boxScatterFitPlotWrap(vut, 'rf_ecc', 'dPrime', ...
    false, false, 6);