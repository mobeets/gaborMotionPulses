%% load data

if figure5_reGenData
    vut = [vuMT_var vu(~[vu.isCell])];
%     vut = vu;
    [scs, scsP, ~] = tools.decodeWithCells(vut, false, false);
    [~, ~, scsA] = tools.decodeWithCells(vut, true, false);
    scsP2 = tools.decodeWithCellsAndShuffle(scsP, 10);
%     save('data/decodeScs3.mat', 'scs', 'scsP');
%     save('data/allCells3.mat', 'scsA');
%     save('data/decodeScsShuffled.mat', 'scsP2');    

else
    x = load('data/allCells3.mat');
    scsA = x.scsA;
    x = load('data/decodeScs3.mat');
    scsP = x.scsP;
    scsP = scsP(~strcmp({scsP.dt}, '20150304a'));
end

Z = num2cell((scsP2(:,1) - mean(scsP2(:,2:end),2)));
[scsP.scoreGainWithShuffle] = Z{:};

%% cellScore vs. rfCorr (same-sign and opposite-sign rfCorr)

Xnm = 'rfDist_norm';
% Xnm = 'noiseCorrLow';
% Xnm = 'noiseCorrAR';
% Ynm = 'pairImprovement';
Ynm = 'pairImprovePctLeft';
ys = [scsP.(Ynm)];
miny = min(ys)-0.01;
maxy = max(ys)+0.01;

ys = [scsP.(Ynm)];
ix = abs(zscore([scsP.(Ynm)])) < 2.5;
scsP0 = scsP(ix);

% ixc = [scsP0.cell1_sep] > 0.6 & [scsP0.cell2_sep] > 0.6;
sum(~ixc)
scsP0 = scsP0(ixc);

% scsP0 = scsP([scsP.rfDist] <= 1);
% ix = sign([scsP0.rfCorr])<0;
ix = [scsP0.cell1_targPref]==[scsP0.cell2_targPref];

% plot.boxScatterFitPlotWrap(scsP0, Xnm, Ynm, true);
% xlabel(Xnm);
% ylim([miny maxy]);
 
plot.boxScatterFitPlotWrap(scsP0(ix), Xnm, Ynm, true);
xlabel([Xnm ' (same sign rfCorr)']);
ylim([miny maxy]);

plot.boxScatterFitPlotWrap(scsP0(~ix), Xnm, Ynm, true);
xlabel([Xnm ' (diff sign rfCorr)']);
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
