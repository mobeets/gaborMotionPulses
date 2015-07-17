%% load data

if figure5_reGenData
    [~, scsP, ~] = tools.decodeWithCells(vu, false, false);
    [~, ~, scsA] = tools.decodeWithCells(vu, true, false);
    save('data/decodeScs2.mat', 'scs', 'scsP');
    save('data/allCells2.mat', 'scsA');
else
    x = load('data/allCells2.mat');
    scsA = x.scsA;
    x = load('data/decodeScs2.mat');
    scsP = x.scsP;
    scsP = scsP(~strcmp({scsP.dt}, '20150304a'));
end

% remove outliers in improvement
Ynm = 'pairImprovement';
ys = [scsP.(Ynm)];
scsP = scsP(~isnan(ys));
ix = zscore([scsP.(Ynm)]) < 3;
scsP = scsP(ix);

%% cellScore vs. rfCorr (same-sign and opposite-sign rfCorr)

Xnm = 'rfDist';
Ynm = 'pairImprovement';
ys = [scsP.(Ynm)];
miny = min(ys)-0.01;
maxy = max(ys)+0.01;

plot.boxScatterFitPlotWrap(scsP(sign([scsP.rfCorr])<0), Xnm, Ynm, false);
xlabel('rfDist (opposite sign rfCorr)');
ylim([miny maxy]);

plot.boxScatterFitPlotWrap(scsP(sign([scsP.rfCorr])>0), Xnm, Ynm, false);
xlabel('rfDist (same sign rfCorr)');
ylim([miny maxy]);

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
