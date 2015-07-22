Ynm = 'pairImprovement';
% Ynm = 'pairImprovePctLeft';
scsP0 = scsP;%(abs(zscore([scsP.(Ynm)]))<2);
Y = [scsP0.(Ynm)];

% scsP0 = scsP;
% Y = -Z;

x0 = [scsP0.singleScoreMax];
x1 = [scsP0.dPrimeMean];

x2 = [scsP0.rfDist];
x2a = [scsP0.rfDist2];
x3 = [scsP0.rfCorr];

x4 = [scsP0.cell1_targPref]==[scsP0.cell2_targPref];
x5 = sign([scsP0.rfCorr])>0;

x6 = [scsP0.noiseCorrAR];
x7 = [scsP0.mnkScore];
x8 = [scsP0.scoreGainWithShuffle];

X = [x0; x1; x8];
% X = [x1; x4; x3];
% X = [x3; x5; x7];

m = fitlm(X', Y')
% figure; plot(m.predict(X')); hold on; plot(Y, 'k');

%%

vs = scsP(zscore(x6) < -2);
for ii = 1:numel(vs)
    figure; hold on;
%     ix = vs(ii).stim == 1;
    v0 = vuMT(strcmp({vuMT.name}, vs(ii).cell1));
    v1 = vuMT(strcmp({vuMT.name}, vs(ii).cell2));
    ix = v0.dirstrength > 0;
    
    scatter(v0.Y(ix), v1.Y(ix), 'r', 'filled');
    scatter(v0.Y(~ix), v1.Y(~ix), 'b', 'filled');
%     scatter(vs(ii).cell1_Y(ix), vs(ii).cell2_Y(ix), 'r', 'filled');
%     scatter(vs(ii).cell1_Y(~ix), vs(ii).cell2_Y(~ix), 'b', 'filled');
    
end

%%

muFcn = @(v) v.wfSvd_U(:,1)*sign(sum(v.wfSvd_U(:,1)));

for ii = 1:numel(vutVar)
    v = vutVar(ii);
    plot.plotKernel2(muFcn(v), zscore(v.Xxy));
    hold on;
    title(v.name);
    scatter(v.rf_center(:,1), v.rf_center(:,2), 300, 'r', 'LineWidth', 2);
    scatter(v.rf_center2(:,1), v.rf_center2(:,2), 300, 'g', 'LineWidth', 2);
    set(gcf, 'Position', [0 0 600 600])
    plot.saveFigure([v.dt '-' v.label], 'rfCenters_MT');
    if mod(ii, 10) == 0
        close all;
    end
end

%%

for ii = 1:numel(scsP)
    v0 = vut(strcmp({vut.name}, scsP(ii).cell1));
    v1 = vut(strcmp({vut.name}, scsP(ii).cell2));
    
    ncA = v0.rf_center;
    ncB = v1.rf_center;
    scsP(ii).rfDist3 = sqrt((ncB(1)-ncA(1))^2 + ...
                    (ncB(2) - ncA(2))^2);
    Ynm = 'YresAR';
    minTrials = 30;
    
    Y1 = v0.(Ynm);
    Y2 = v1.(Ynm);
    ix = ~isnan(Y1) & ~isnan(Y2);
    if sum(ix) < minTrials
        nc = nan;
        return;
    end
    Y1 = Y1(ix);
    Y2 = Y2(ix);
    nc = corr(Y1, Y2);

    scsP(ii).noiseCorrAR = corr(Y1, Y2);

end
