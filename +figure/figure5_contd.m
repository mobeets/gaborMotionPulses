
Ynm = 'scoreGainWithShuffle';

ix = abs(zscore([scsP.(Ynm)])) < 2;
figure; hist(zscore([scsP.(Ynm)]));
scsP0 = scsP(ix);

%%

ix = [scsP0.sameTarg];
ix1 = [scsP0.sameCorr];
ix2a = (ix&~ix1); % same dir, anti-corr
ix2b = (~ix&ix1); % opp dir, pos-cor
ix2 = ix2a | ix2b;

xnm = 'rfDist3';
% xnm = 'rfCorr';
% xnm = 'noiseCorrAR';
% ynm = 'noiseCorrAR';
ynm = 'rfCorr';
% ynm = 'scoreGainWithShuffle';

figure; hold on; set(gca, 'FontSize', 14); set(gcf, 'color', 'w');
t1 = scsP0(ix2a);
t2 = scsP0(ix2b);
t3 = scsP0(~ix2 & ~ix);
t4 = scsP0(~ix2 & ix);
scatter([t3.(xnm)], [t3.(ynm)], 80, 'r');
scatter([t1.(xnm)], [t1.(ynm)], 80, 'b');
scatter([t2.(xnm)], [t2.(ynm)], 'r', 'filled');
scatter([t4.(xnm)], [t4.(ynm)], 'b', 'filled');
h = legend('diff       dir, (-) sig-corr', 'same dir, (-) sig-corr', 'diff       dir, (+) sig-corr', 'same dir, (+) sig-corr');
set(h, 'Location', 'NorthEastOutside');
set(gcf, 'Position', [600 600 800 400]);
xlabel(xnm);
ylabel(ynm);

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