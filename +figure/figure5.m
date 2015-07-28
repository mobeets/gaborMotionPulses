%% load data

if figure5_reGenData
    vut = [vuMT_var vu(~[vu.isCell])];
%     vut = vuMT;
    [scs, scsP, ~] = tools.decodeWithCells(vut, false, false);
    [scs0, scs1] = tools.decodeWithCellsAndShuffle({scsP.stim}, ...
        {scsP.Ys}, 10);
    [~, ~, scsA] = tools.decodeWithCells(vut, true, false);    
    [scs3, scs4] = tools.decodeWithCellsAndShuffle({scsA.stim}, ...
        {scsA.Ys}, 10);
    save('data/decodeScs_all.mat', 'scs', 'scsP');
    save('data/allCells_all.mat', 'scsA');
    save('data/decodeScsShuffled_scsP_all.mat', 'scs0', 'scs1');
    save('data/decodeScsShuffled_scsA_all.mat', 'scs3', 'scs4');

else
    x = load('data/allCells.mat');
    scsA = x.scsA;
    x = load('data/decodeScs.mat');
    scsP = x.scsP;
    x = load('data/decodeScsShuffled_scsP.mat');
    scs0 = x.scs0; scs1 = x.scs1;
    x = load('data/decodeScsShuffled_scsA.mat');
    scs3 = x.scs3; scs4 = x.scs4;
end

%%

vut = [vuMT_var vu(~[vu.isCell])];
[scs, scsP, ~] = tools.decodeWithCells(vut, false, false);
[scs0, scs1] = tools.decodeWithCellsAndShuffle({scsP.stim}, ...
    {scsP.Ys}, 10);
save('data/decodeScs_all_2.mat', 'scs', 'scsP');
save('data/decodeScsShuffled_scsP_all_2.mat', 'scs0', 'scs1');

%% sorry this is messy...

Z1a = scs0;
Z1b = mean(scs1, 2);
Z1 = Z1a - Z1b;
Z1s = (std(scs1')/sqrt(size(scs1,2)))';
Z1_lb = Z1a - (Z1b+2*Z1s);
Z1_ub = Z1a - (Z1b-2*Z1s);

Z1 = num2cell(Z1);
Z1a = num2cell(Z1a);
Z1b = num2cell(Z1b);
Z1s = num2cell(Z1s);
Z1_lb = num2cell(Z1_lb);
Z1_ub = num2cell(Z1_ub);

[scsP.scoreGainWithCorrs] = Z1{:};
[scsP.scoreNoShuffle] = Z1a{:};
[scsP.scoreWithShuffle] = Z1b{:};
[scsP.scoreGainWithCorrs_lb] = Z1_lb{:};
[scsP.scoreGainWithCorrs_ub] = Z1_ub{:};

Z1a = scs0;
Z1b = mean(scs1, 2);
Z1 = Z1a - Z1b;
Z1s = (std(scs1')/sqrt(size(scs1,2)))';
Z1_lb = Z1a - (Z1b+2*Z1s);
Z1_ub = Z1a - (Z1b-2*Z1s);

Z1 = num2cell(Z1);
Z1a = num2cell(Z1a);
Z1b = num2cell(Z1b);
Z1s = num2cell(Z1s);
Z1_lb = num2cell(Z1_lb);
Z1_ub = num2cell(Z1_ub);

[scsA.scoreGainWithCorrs] = Z1{:};
[scsA.scoreNoShuffle] = Z1a{:};
[scsA.scoreWithShuffle] = Z1b{:};
[scsA.scoreGainWithCorrs_lb] = Z1_lb{:};
[scsA.scoreGainWithCorrs_ub] = Z1_ub{:};

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

ixc = [scsP0.cell1_sep] > 0.6 & [scsP0.cell2_sep] > 0.6;
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
figure; hold on; set(gca, 'FontSize', 20); set(gcf, 'color', 'white');
for ii = 1:numel(scsA)
    cix = scsA(ii).ncells;
    scatter(scsA(ii).mnkScore, scsA(ii).score, 100, clrs(cix,:), 'filled');
    scatter(scsA(ii).mnkScore, scsA(ii).score, 100, 'k');
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
