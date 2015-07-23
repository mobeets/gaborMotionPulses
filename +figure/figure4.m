
Yh = 'YresAR';
cpY = 'cp_YresAR';

% vut = vuMT([vuMT.nzertrials] >= 30);
% vut = vuMT(arrayfun(@(v) all(v.nlowmottrials > 30), vuMT));
vut = vuMT_var;

%% noise-correlation vs. distance of RF centers
% plot.pairwiseCorrs(vut, 'rf_center', Yh, ...
%     nan,nan,nan,@(x,y) norm(x-y,2));
% set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
% plot.saveFigure('MT - noise-corr vs. rf-center', figDir, gcf, 'pdf');

vut2 = vut([vut.separability_index]>0.6);
% getWf = @(v) v.wfSvd_U(:,1)*sign(sum(v.wfSvd_U(:,1)));
getWf = @(v) v.w(:);
sameTarg = @(v0, v1) v0.targPref == v1.targPref;
sameCorr = @(v0, v1) corr(getWf(v0), getWf(v1)) > 0;
[pss, mdl, nms] = plot.pairwiseCorrs(vut2, 'rf_center', Yh, ...
    nan, 'dist(rfCenter_1, rfCenter_2) for same targPref', ...
    nan, @(x,y) norm(x-y,2), @corr, sameCorr);

diffTarg = @(v0, v1) v0.targPref ~= v1.targPref;
diffCorr = @(v0, v1) corr(getWf(v0), getWf(v1)) < 0;
[pss, mdl, nms] = plot.pairwiseCorrs(vut2, 'rf_center', Yh, ...
    nan, 'dist(rfCenter_1, rfCenter_2) for opposite targPref', ...
    nan, @(x,y) norm(x-y,2), @corr, diffCorr);

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

%%

plot.pairwiseCorrs(vut, 'w', 'YresAR', ...
    nan, 'jsDivergence', ...
    nan, @(v0,v1) tools.jsDivergence(v0, v1, true), @corr);

%%

plot.pairwiseCorrs(vut, 'wrat', 'YresAR', ...
    nan, 'wrat', ...
    nan, @(a,b) a+b, @corr, diffTarg);
plot.pairwiseCorrs(vut, 'wrat', 'YresAR', ...
    nan, 'wrat', ...
    nan, @(a,b) a+b, @corr, sameTarg);
