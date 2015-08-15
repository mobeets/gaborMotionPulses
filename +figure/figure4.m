
Yh = 'YresAR';
cpY = 'cp_YresAR';

pyS = struct();

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
[a,b,c] = plot.pairwiseCorrs(vut, 'w', Yh);
set(gcf, 'PaperSize', [4 6], 'PaperPosition', [0 0 4 6])
plot.saveFigure('MT - noise-corr vs. rf-corr', figDir, gcf, 'png');

%% CP vs. dPrime

m = plot.boxScatterFitPlotWrap(vut, 'dPrime', cpY);
plot.saveFigure(['MT - ' cpY ' vs. dPrime'], figDir, gcf);

%% RF eccentricity vs. CP
plot.boxScatterFitPlotWrap(vut, 'rf_ecc', cpY, ...
    true, false, 6);
set(gca, 'Xtick', [0 .5 1 1.5])
hold on
plot([0 1.5], [.5 .5], 'k:')
figure.cleanupForPrint(gcf, 'PaperSize', [25 25], 'FontSize', 8)
plot.saveFigure(['MT - ' cpY ' vs. rf-eccentricity'], figDir, gcf);

%% output struct
pyS.rf_ecc = [vut.rf_ecc];
pyS.rf_xy  = reshape([vut.rf_center], 2,[]);
pyS.dprime = [vut.dPrime];
pyS.CPres  = [vut.cp_YresAR];
pyS.CPlow  = [vut.cp_Ylow];

save(fullfile(figDir, 'figure04_pystruct.mat'), '-v7.3', '-struct', 'pyS')
    
%% dPrime vs. rfEcc
plot.boxScatterFitPlotWrap(vut, 'rf_ecc', 'dPrime', ...
    true, false, 6);
set(gca, 'Xtick', [0 .5 1 1.5])
hold on
figure.cleanupForPrint(gcf, 'PaperSize', [25 25], 'FontSize', 8)
axis off
title('')
plot.saveFigure(['MT - ' 'dPrime' ' vs. rf-eccentricity'], figDir, gcf);

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

%% Figure 4 Jake Draft
showBox=false;
showPCline=true;
figure(4); clf
subplot(2,2,1)
[~,xs,ys]=plot.boxScatterFitPlotWrap(vut, 'rf_ecc', 'dPrime', showBox, false, 10, showPCline);
[~,idx] = sort(xs);
plot(xs(idx), smooth(xs(idx), ys(idx), 0.3, 'loess'), 'r-');
title('Sensitivity')
xlabel('RF Center Distance')
ylabel('d''')
git pull
subplot(2,2,2)
[~,xs,ys]=plot.boxScatterFitPlotWrap(vut, 'rf_ecc', cpY, showBox, false, 10, showPCline);
[~,idx] = sort(xs);
plot(xs(idx), smooth(xs(idx), ys(idx), 0.3, 'loess'), 'r-');
xlim([0 1.2])
ylim([.4 .65])
title('Choice Probability')
xlabel('RF Center Distance')
ylabel('CP')
plot(xlim, [.5 .5], 'k')




