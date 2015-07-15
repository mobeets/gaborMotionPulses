vv = vw(~[vw.isLIP] & ([vw.dPrime] > 0.2 | ~[vw.isCell]));
[scs, scsP, scsA] = decode.main(vv);

%%

vv = vw(~[vw.isMT] & ([vw.dPrime] > 0.08 | ~[vw.isCell]));
[scs, scsP, scsA] = decode.main(vv);

%%

clrs = cbrewer('seq', 'Blues', max([scsA.ncells]), 'pchip');
figure; hold on; set(gca, 'FontSize', 14); set(gca, 'color', 'white');
for ii = 1:numel(scsA)
    cix = scsA(ii).ncells;
    scatter(scsA(ii).mnkScore, scsA(ii).score, 50, clrs(cix,:), 'filled');
    scatter(scsA(ii).mnkScore, scsA(ii).score, 50, 'k');
end
xlabel('mnkScore');
ylabel('score');
xlim(ylim);
plot(xlim, xlim, 'k--');

figure; hold on; set(gca, 'FontSize', 14); set(gca, 'color', 'white');
for ii = 1:numel(scsA)
    cix = scsA(ii).ncells;
    scatter(scsA(ii).ncells, scsA(ii).score-scsA(ii).mnkScore, 50, clrs(cix,:), 'filled');
    scatter(scsA(ii).ncells, scsA(ii).score-scsA(ii).mnkScore, 50, 'k');
end
xlabel('ncells');
ylabel('score');

%%

scsP = scsP(~strcmp({scsP.dt}, '20150304a'));
scs = scs(~strcmp({scs.dt}, '20150304a'));

%%

Ynm = 'pairImprovement';
Ynm = 'pairImprovePctLeft';
% Ynm = 'score';

%%

plot.boxScatterFitPlotWrap(scsP, 'singleScoreMax', 'score', false, false)
xlim(ylim); plot(xlim, xlim, 'k--');
plot.boxScatterFitPlotWrap(scsP, 'mnkScore', 'score', false, false)
xlim(ylim); plot(xlim, xlim, 'k--');


%%

plot.boxScatterFitPlotWrap(scsP, 'singleScoreMax', Ynm)
plot.boxScatterFitPlotWrap(scsP, 'dPrimeMax', Ynm)

plot.boxScatterFitPlotWrap(scsP, 'dPrimeMax', 'singleScoreMax')

%%

plot.boxScatterFitPlotWrap(scsP, 'rfCorr', Ynm, false, false)
%%
plot.boxScatterFitPlotWrap(scsP, 'noiseCorrZer', Ynm)
plot.boxScatterFitPlotWrap(scsP, 'noiseCorrLow', Ynm)
%%
plot.boxScatterFitPlotWrap(scsP, 'rfCorr', 'noiseCorrZer')
plot.boxScatterFitPlotWrap(scsP, 'rfCorr', 'noiseCorrLow')

%%

plot.boxScatterFitPlotWrap(scsP, 'noiseCorrZerRfCorr', Ynm)
plot.boxScatterFitPlotWrap(scsP, 'noiseCorrLowRfCorr', Ynm)
plot.boxScatterFitPlotWrap(scsP, 'noiseCorrZerRfCorr2', Ynm)
plot.boxScatterFitPlotWrap(scsP, 'noiseCorrLowRfCorr2', Ynm)


%%

Ynm = 'pairImprovement';
miny = min([scsP.(Ynm)])-0.01;
maxy = max([scsP.(Ynm)])+0.01;

ys = [scsP.(Ynm)];
S = scsP(~isnan(ys));
ix = zscore([S.(Ynm)]) < 3;

plot.boxScatterFitPlotWrap(S(ix & sign([S.rfCorr])<0), 'rfDist', Ynm, false)
ylim([miny maxy]);
plot.boxScatterFitPlotWrap(S(ix & sign([S.rfCorr])>0), 'rfDist', Ynm, false)
ylim([miny maxy]);

%%

plot.boxScatterFitPlotWrap(scsP, 'rfDist', 'noiseCorrZer')
plot.boxScatterFitPlotWrap(scsP, 'rfDist', 'noiseCorrLow')

