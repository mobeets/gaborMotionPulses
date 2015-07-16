%% load fit data

fitdir = '20150615';
fitbasedir = 'data';
ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
vu = [vp vn];

%% pairwise plots

Yh = 'Yres';
% vut = vu([vu.nzertrials] >= 30);

% noise-correlation vs. distance of RF centers
plot.pairwiseCorrs(vu([vu.isMT]), 'rf_center', Yh, ...
    nan,nan,nan,@(x,y) norm(x-y,2));
% plot.pairwiseCorrs(vut([vut.isMT]), 'rf_center', Yh, ...
%     nan,nan,nan,@(x,y) norm(x-y,2));
plot.saveFigure('MT - noise-corr vs. rf-center', 'tmp', gcf);



% noise-correlation vs. RF similarity
% plot.pairwiseCorrs(vu([vu.isLIP]), 'w', Yh);
% plot.saveFigure('LIP - noise-corr vs. rf-corr', 'tmp', gcf);
plot.pairwiseCorrs(vu([vu.isMT]), 'w', Yh);
% plot.pairwiseCorrs(vut([vut.isMT]), 'w', Yh);
plot.saveFigure('MT - noise-corr vs. rf-corr', 'tmp', gcf);

%% CP plots

cpY = 'cp_YresARc';
% cpY = 'cp_Ylow'
% vut = vu;
% vut = vu([vu.nzertrials] >= 30);
% vut = vu(arrayfun(@(v) all(v.nlowmottrials > 30), vu));
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
