
Ynm = 'scoreGainWithCorrs';

ix = abs(zscore([scsP.(Ynm)])) < 2;
% figure; hist(zscore([scsP.(Ynm)]));
scsP0 = scsP(ix);

scsP0 = scsP;

% 
% scsP0 = scsP0([scsP0.noiseCorr_pavg] < 0.05);

%%

ix3 = [scsP0.sameTarg] & ... % same pool
    [scsP0.noiseCorrAR] > 0 & ... % positive noise corr
    [scsP0.scoreGainWithCorrs_lb] > 0; % but better than shuffle?
ix4 = ~[scsP0.sameTarg] & ... % opposite pools
    [scsP0.noiseCorrAR] < 0 & ...% negative noise corr
    [scsP0.scoreGainWithCorrs_lb] > 0; % but better than shuffle?
t1 = scsP0(ix3);
t2 = scsP0(ix4);

for ii = 1:numel(t2)
%     plot.cellPairKernelsSingle(vuMT, t2(ii).cell1, t2(ii).cell2);
    v = t2(ii);
    figure; S = plot.visualizePairwiseCorr(v.Ys, v.stim);
    title([v.cell1 ', ' v.cell2 ': ' num2str(S.pcOpt - S.pcBlind)]);
end
%%

figure; hold on; set(gca, 'FontSize', 14); set(gcf, 'color', 'w');
% xnm = 'rfDist_norm';
% xnm = 'rfCorr';
xnm = 'noiseCorrAR';
% ynm = 'noiseCorrAR';
% ynm = 'rfCorr';
ynm = 'scoreGainWithCorrs_lb';
% ynm = 'scoreOverDiagLinear';
znm = 'rfCorr';

c1 = [0.9 0.7 0.2];
c2 = [0.3 0.6 0.6];
scatter3([t1.(xnm)], [t1.(ynm)], [t1.(znm)], 50, c1, 'filled');
scatter3([t2.(xnm)], [t2.(ynm)], [t2.(znm)], 50, c2, 'filled');
xlabel(xnm);
ylabel(ynm);
zlabel(znm);
plot(xlim, [0 0], 'k--');
plot([0 0], ylim, 'k--');

%%
ix = [scsP0.sameTarg];
ix1 = [scsP0.sameCorr];
ix2a = (ix&~ix1); % same dir, anti-corr
ix2b = (~ix&ix1); % opp dir, pos-cor
ix2 = ix2a | ix2b;

% xnm = 'rfDist_norm';
% xnm = 'rfCorr';
% xnm = 'sigNoiseAngleDev';
% xnm = 'noiseCorr_slope';
% xnm = 'noiseCorr_avg';
xnm = 'noiseCorrAR';

% ynm = 'noiseCorrAR';
% ynm = 'noiseCorr_avg';
% ynm = 'rfCorr';
ynm = 'scoreGainWithCorrs_lb';
% ynm = 'scoreOverDiagLinear';
znm = 'rfCorr';

figure; hold on; set(gca, 'FontSize', 14); set(gcf, 'color', 'w');
t1 = scsP0(ix2a); % same dir, anti-corr
t2 = scsP0(ix2b); % opp dir, pos-cor
t3 = scsP0(~ix2 & ~ix);
t4 = scsP0(~ix2 & ix);
c2 = [0.3 0.6 0.6];
c1 = [0.9 0.7 0.2];
scatter3([t3.(xnm)], [t3.(ynm)], [t3.(znm)], 50, c1, 'filled');
scatter3([t2.(xnm)], [t2.(ynm)], [t2.(znm)], 80, c1);
scatter3([t1.(xnm)], [t1.(ynm)], [t1.(znm)], 80, c2);
scatter3([t4.(xnm)], [t4.(ynm)], [t4.(znm)], 50, c2, 'filled');
h = legend('diff dir, (-) corr', 'diff dir, (+) corr', 'same dir, (-) corr', 'same dir, (+) corr');
set(h, 'Location', 'NortheastOutside');
set(gcf, 'Position', [600 600 800 400]);
xlabel(xnm);
ylabel(ynm);
zlabel(znm);
plot(xlim, [0 0], 'k--');
plot([0 0], ylim, 'k--');

%%

Ynm = 'noiseCorrAR';
figure; hold on; set(gca, 'FontSize', 14); set(gcf, 'color', 'w');
bins = linspace(-1, 1, 11);
subplot(2,2,1); hold on;
hist([scsP0(ix2a).(Ynm)], bins);
title('same dir, (-) sig-corr');
subplot(2,2,2); hold on;
hist([scsP0(ix2b).(Ynm)], bins);
title('opp dir, (+) sig-corr');
subplot(2,2,3); hold on;
hist([scsP0(~ix2 & ~ix).(Ynm)], bins);
title('opp dir, (-) sig-corr');
subplot(2,2,4); hold on;
hist([scsP0(~ix2 & ix).(Ynm)], bins);
title('same dir, (+) sig-corr');

%%
tm = scsP0(ix2);
plot.boxScatterFitPlotWrap(tm, 'noiseCorrAR', Ynm);

%%

plot.boxScatterFitPlotWrap(scsP0(ix), 'noiseCorrAR', Ynm);
xlabel('noise-corr for same targPref');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'noiseCorrAR', Ynm);
xlabel('noise-corr for diff targPref');
%%
plot.boxScatterFitPlotWrap(scsP0, 'rfDist3', Ynm);
plot.boxScatterFitPlotWrap(scsP0(ix), 'rfDist', Ynm);
xlabel('dist(rfCenter_1, rfCenter_2) for same targPref');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfDist3', Ynm);
xlabel('dist(rfCenter_1, rfCenter_2) for diff targPref');

%%

% plot.boxScatterFitPlotWrap(scsP0, 'rfCorr', 'scoreGainWithShuffle');
plot.boxScatterFitPlotWrap(scsP0(ix), 'rfCorr', 'scoreGainWithShuffle');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfCorr', 'scoreGainWithShuffle');

%%

plot.boxScatterFitPlotWrap(scsP0(ix), 'rfCorr', 'rfDist_norm');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfCorr', 'rfDist_norm');

plot.boxScatterFitPlotWrap(scsP0, 'rfDist_norm', 'scoreGainWithShuffle');
% xlabel('dist(rfCenter_1, rfCenter_2) for same targPref');
% plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfCorr', Ynm);
% xlabel('dist(rfCenter_1, rfCenter_2) for diff targPref');

%%
plot.boxScatterFitPlotWrap(scsP0(ix), 'rfDist_norm', 'scoreGainWithShuffle');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfDist_norm', 'scoreGainWithShuffle');
