%% load data

if figure5_reGenData
    vut = [vuMT_var vu(~[vu.isCell])];
%     vut = vu;
    [scs, scsP, ~] = tools.decodeWithCells(vut, false, false);
    [~, ~, scsA] = tools.decodeWithCells(vut, true, false);
%     save('data/decodeScs3.mat', 'scs', 'scsP');
%     save('data/allCells3.mat', 'scsA');
else
    x = load('data/allCells3.mat');
    scsA = x.scsA;
    x = load('data/decodeScs3.mat');
    scsP = x.scsP;
    scsP = scsP(~strcmp({scsP.dt}, '20150304a'));
end

% remove outliers in improvement
% Ynm = 'pairImprovement';
% ys = [scsP.(Ynm)];
% scsP = scsP(~isnan(ys));
% ix = zscore([scsP.(Ynm)]) < 3;
% scsP = scsP(ix);

%% cellScore vs. rfCorr (same-sign and opposite-sign rfCorr)

Xnm = 'rfDist3';
% Xnm = 'noiseCorrLow';
Ynm = 'pairImprovement';
% Ynm = 'pairImprovePctLeft';
ys = [scsP.(Ynm)];
miny = min(ys)-0.01;
maxy = max(ys)+0.01;

ys = [scsP.(Ynm)];
ix = abs(zscore([scsP.(Ynm)])) < 2;
scsP0 = scsP(ix);

% scsP0 = scsP([scsP.rfDist] <= 1);
% ix = sign([scsP0.rfCorr])<0;
% ix = [scsP.cell1_targPref]==[scsP.cell2_targPref];

plot.boxScatterFitPlotWrap(scsP0, Xnm, Ynm, true);
xlabel(Xnm);
ylim([miny maxy]);
% 
% plot.boxScatterFitPlotWrap(scsP(ix), Xnm, Ynm, true);
% xlabel([Xnm ' (opposite sign rfCorr)']);
% ylim([miny maxy]);
% 
% plot.boxScatterFitPlotWrap(scsP(~ix), Xnm, Ynm, true);
% xlabel([Xnm ' (same sign rfCorr)']);
% ylim([miny maxy]);

%% cellScore vs. mnkScore

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

%% cellScore vs. nCells

figure; hold on; set(gca, 'FontSize', 14); set(gca, 'color', 'white');
for ii = 1:numel(scsA)
    cix = scsA(ii).ncells;
    scatter(scsA(ii).ncells, scsA(ii).score, 50, clrs(cix,:), 'filled');
    scatter(scsA(ii).ncells, scsA(ii).score, 50, 'k');
end
xlabel('ncells');
ylabel('score');
