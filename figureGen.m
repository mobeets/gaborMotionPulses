%% load fit data

fitdir = '20150615';
fitbasedir = 'data';
ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
vu = [vp vn];

%% pairwise plots

% noise-correlation vs. distance of RF centers
plot.pairwiseCorrs(vu([vu.isMT]), 'rf_center', 'Yfrz', ...
    nan,nan,nan,@(x,y) norm(x-y,2));
plot.saveFigure('MT - noise-corr vs. rf-center', 'tmp', gcf);

% noise-correlation vs. RF similarity
plot.pairwiseCorrs(vu([vu.isLIP]), 'w', 'Yfrz');
plot.saveFigure('LIP - noise-corr vs. rf-corr', 'tmp', gcf);
plot.pairwiseCorrs(vu([vu.isMT]), 'w', 'Yfrz');
plot.saveFigure('MT - noise-corr vs. rf-corr', 'tmp', gcf);

%% CP plots

vut = vu([vu.nzertrials] >= 30);

plot.boxScatterFitPlotWrap(vut([vut.isMT]), 'dPrime', 'cp_Yzer');
plot.saveFigure('MT - cp_Yzer vs. dPrime', 'tmp', gcf);

plot.boxScatterFitPlotWrap(vut([vut.isMT]), 'rf_ecc', 'cp_Yzer', ...
    true, true, 6);
xlim([0 0.6]);
plot.saveFigure('MT - cp_Yzer vs. rf-eccentricity', 'tmp', gcf);

plot.boxScatterFitPlotWrap(vu([vu.isMT]), 'dPrime', 'rf_ecc', ...
    false, false, 6);
