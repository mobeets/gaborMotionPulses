scs = shuffleDecode(scsP, 20);
sc = cell2mat(scs);
%%
Z = (sc(:,1) - mean(sc(:,2:end),2));
Z2 = Z./(1-sc(:,1));
for ii = 1:numel(scsP)
    scsP(ii).scoreGainWithShuffle = Z(ii);
    scsP(ii).scoreGainWithShuffle2 = Z2(ii);
    scsP(ii).sameTarg = scsP(ii).cell1_targPref == scsP(ii).cell2_targPref;
    scsP(ii).sameCorr = sign(scsP(ii).rfCorr)>0;
end
%%

Ynm = 'scoreGainWithShuffle';

ix = abs(zscore([scsP.(Ynm)])) < 2;
figure; hist(zscore([scsP.(Ynm)]));
scsP0 = scsP(ix);

ix = [scsP0.sameTarg];

plot.boxScatterFitPlotWrap(scsP0(ix), 'noiseCorrAR', Ynm);
xlabel('noise-corr for same targPref');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'noiseCorrAR', Ynm);
xlabel('noise-corr for diff targPref');

plot.boxScatterFitPlotWrap(scsP0, 'rfDist', Ynm);
plot.boxScatterFitPlotWrap(scsP0(ix), 'rfDist', Ynm);
xlabel('dist(rfCenter_1, rfCenter_2) for same targPref');
plot.boxScatterFitPlotWrap(scsP0(~ix), 'rfDist', Ynm);
xlabel('dist(rfCenter_1, rfCenter_2) for diff targPref');
