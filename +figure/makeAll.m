%% init

fitnm = 'space-time';
doSaveFigs = false;

saveDir = fullfile('data', 'figures', fitnm);
if doSaveFigs && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
pairsFnm = fullfile(saveDir, 'allPairs.mat');
decodingFnm = fullfile(saveDir, 'allPairsWithScs.mat');
% addpath(genpath('~/code/gaborMotionASD/mASD/'));
% or: mpm install masd -u https://github.com/mobeets/masd.git

%% make fits

% fits
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');
% fitAllSTRFs('time-only', false, 'cells ASD time-only', io.getDates, 'MT');
% fitAllSTRFs('flat', false, 'cells Flat', io.getDates, 'MT');

%% decode with shuffled cell pairs, triplets, quads, etc.

% decoding analysis, with noise corr shuffle
cellGroups = tools.makeCellGroups(cells, [2 3 4]);

% nShuffles = 1;
% scs = tools.decodeAndShuffle({cellGroups.stimdir}, ...
%     {cellGroups.Ys}, nShuffles);
% scDelta = num2cell(scs.scsDelta);
% [cellGroups.scoreGainWithCorrs] = scDelta{:};
% scDeltaSe = num2cell(scs.scsShufSe);
% [cellGroups.scoreGainWithCorrs_se] = scDeltaSe{:};
% 
% save(decodingFnm, 'scs', 'cellGroups');

%% make cell pairs, and perform decoding analysis

% % cell pairs
% allPairs = tools.makeCellPairs(cells);
% save(pairsFnm, 'allPairs');
% 
% % decoding analysis, with noise corr shuffle
% nShuffles = 25;
% scs = tools.decodeAndShuffle({allPairs.stimdir}, ...
%     {allPairs.Ys}, nShuffles);
% scDelta = num2cell(scs.scsDelta);
% [allPairs.scoreGainWithCorrs] = scDelta{:};
% scDeltaSe = num2cell(scs.scsShufSe);
% [allPairs.scoreGainWithCorrs_se] = scDeltaSe{:};
% save(decodingFnm, 'scs', 'allPairs');

%% load cells and cell pairs

cells = tools.makeFitSummaries(['data/fits/' fitnm]);
d = load(decodingFnm, 'scs', 'allPairs');
allPairs = d.allPairs;

ixCellsToKeep = figure.filterCellsAndPairs(cells, false);
ixPairsToKeep = figure.filterCellsAndPairs(allPairs, true);
pairs = allPairs(ixPairsToKeep);

%% summarize experiment, cell, and pair counts

% get all experiments
curCells = cells(ixCellsToKeep);
dts = {curCells.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
sessions = unique(dts);
ncells_per_session = grpstats([curCells.index], dts, @numel);
ntrials_per_session = grpstats([curCells.ntrials], dts, @mean);

disp('-----SUMMARY-----');
% summarize experiments, cells, and trials per monkey
isNancy = sessions >= 20150101;
nexps = [sum(~isNancy) sum(isNancy)];
fprintf('# sessions: P=%d, N=%d\n', nexps);
ncells = [sum(ncells_per_session(~isNancy)) sum(ncells_per_session(isNancy))];
fprintf('# cells: P=%d, N=%d\n', ncells);

mu = [mean(ncells_per_session(~isNancy)) mean(ncells_per_session(isNancy))];
sd = [std(ncells_per_session(~isNancy)) std(ncells_per_session(isNancy))];
avgcells = [mu; sd]; avgcells = avgcells(:);
fprintf('avg # cells per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f\n', avgcells);

mu = [mean(ntrials_per_session(~isNancy)) mean(ntrials_per_session(isNancy))];
sd = [std(ntrials_per_session(~isNancy)) std(ntrials_per_session(isNancy))];
avgtrials = [mu; sd]; avgtrials = avgtrials(:);
fprintf('avg # trials per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f\n', avgtrials);

% summarize pair filtering as well
ixCellsToKeepInPairs = figure.filterCellsAndPairs(cells, true);
curCells = cells(ixCellsToKeepInPairs);
dts = {curCells.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
isNancy = dts >= 20150101;
ncells = [sum(~isNancy) sum(isNancy)];
fprintf('# cells used in pairs: P=%d, N=%d\n', ncells);

dts = {pairs.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
isNancy = dts >= 20150101;
ncells = [sum(~isNancy) sum(isNancy)];
fprintf('# pairs: P=%d, N=%d\n', ncells);

%% Fig 2 - ASD

% cellName = '20140304-MT_12'; % 300-400
cellName = '20150304b-MT_14'; % 600-700; 60-160
startTrial = 600;
endTrial = 700;

% 2b - stimuli
fig2b = figure.plotStimOrSpikes(cells, cellName, startTrial, endTrial, 'stim');

% 2c - spikes
fig2c = figure.plotStimOrSpikes(cells, cellName, startTrial, endTrial, 'spikes');

% 2d - all RFs
fig2d = figure.showAllRfs(cells(ixCellsToKeep));

if doSaveFigs
    curSaveDir = fullfile(saveDir, 'Fig2');
    if ~exist(curSaveDir, 'dir')
        mkdir(curSaveDir);
    end
    export_fig(fig2b, fullfile(curSaveDir, '2b.pdf'));
    export_fig(fig2c, fullfile(curSaveDir, '2c.pdf'));
    export_fig(fig2d, fullfile(curSaveDir, '2d.pdf'));
end

%% Fig 3 - decoding

doSaveFigs = false;

% 3b
fig3b = figure.deltaDecodingAcc(pairs);

% 3a, 3c, 3e
xnm = 'rfCorr';
ynm = 'noiseCorrAR';
fig3a = figure.deltaDecodingScatters(pairs, xnm, ynm, 'both');
fig3c = figure.deltaDecodingScatters(pairs, xnm, ynm, 'same');
fig3e = figure.deltaDecodingScatters(pairs, xnm, ynm, 'different');

% 3d, 3f
xnm = 'noiseCorrAR';
ynm = 'scoreGainWithCorrs';
fig3d = figure.deltaDecodingScatters(pairs, xnm, ynm, 'same');
fig3f = figure.deltaDecodingScatters(pairs, xnm, ynm, 'different');

if doSaveFigs
    curSaveDir = fullfile(saveDir, 'Fig3');
    if ~exist(curSaveDir, 'dir')
        mkdir(curSaveDir);
    end
    export_fig(fig3a, fullfile(curSaveDir, '3a.pdf'));
    export_fig(fig3b, fullfile(curSaveDir, '3b.pdf'));
    export_fig(fig3c, fullfile(curSaveDir, '3c.pdf'));
    export_fig(fig3d, fullfile(curSaveDir, '3d.pdf'));
    export_fig(fig3e, fullfile(curSaveDir, '3e.pdf'));
    export_fig(fig3f, fullfile(curSaveDir, '3f.pdf'));
end

%% summarize cell pair types

disp('---------');
scs = pairs;
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
xnm = 'noiseCorrAR'; xnmA = 'r_{sc} >= 0'; xnmB = 'r_{sc} < 0';
ynm = 'scoreGainWithCorrs'; ynmA = '\Delta >= 0'; ynmB = '\Delta < 0';
% xnm = 'sameTarg'; xnmA = 'same pools'; xnmB = 'diff pools';
% ynm = 'sameCorr'; ynmA = 'r_{RF} >= 0'; ynmB = 'r_{RF} < 0';

xG = [scs.(xnm)] > 0;
yG = [scs.(ynm)] > 0;

disp(['# in ' xnmA ' vs. ' xnmB ':']);
[ppct(sum(xG)) ppct(sum(~xG))]
disp(['# in ' ynmA ' vs. ' ynmB ':']);
[ppct(sum(yG)) ppct(sum(~yG))]
disp(['# with ' ynmA ' and ' xnmA ' vs. ' xnmB ':']);
[ppct(sum(xG & yG)) ppct(sum(~xG & yG))]
disp(['# with ' ynmB ' and ' xnmA ' vs. ' xnmB ':']);
[ppct(sum(xG & ~yG)) ppct(sum(~xG & ~yG))]
disp('---------');

% of pairs in same pool, % with r_RF > 0 and r_sc > 0
scs = pairs;
scs = scs([scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] > 0;
B = [scs.noiseCorrAR] > 0;
disp(['# in same pool with r_RF > 0 and r_sc > 0: ' ppct(sum(A & B))]);

% of pairs in diff pools, % with r_RF < 0
scs = pairs;
scs = scs(~[scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] < 0;
% B = [scs.noiseCorrAR] < 0;
disp(['# in diff pools with r_{RF} < 0: ' ppct(sum(A))]);

% of pairs in diff pools, % with r_RF > 0 and r_sc > 0
scs = pairs;
scs = scs(~[scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] > 0;
B = [scs.noiseCorrAR] > 0;
disp(['# in diff pools with r_{RF} > 0 and r_{sc} > 0: ' ppct(sum(A & B))]);

% # with significant \Delta > 0 or \Delta < 0
scs = pairs;
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
S = [scs.scoreGainWithCorrs];
se = [scs.scoreGainWithCorrs_se];
a = (S - se > 0);
b = (S + se < 0);
disp(['# with significant \Delta > 0: ' ppct(sum(a)) ...
    '; \Delta < 0: ' ppct(sum(~b))]);

% correlation between r_sc and r_RF
scs = pairs;
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.noiseCorrAR];
B = [scs.rfCorr];
ix = ~isnan(A) & ~isnan(B);
r = corr(A(ix)', B(ix)');
disp(['correlation between r_{sc} and r_{RF}: ' sprintf('%0.2f', r)]);

% avg noise corr
mu = nanmean([pairs.noiseCorrAR]);
sd = nanstd([pairs.noiseCorrAR]);%/sqrt(numel(pairs));
disp(['avg. noise corr: ' sprintf('%0.2f +/- %0.2f', [mu sd])]);

%% Fig 4 - cell pair examples

figure.plotExamplePairs(pairs, cells, doSaveFigs, ...
    fullfile(saveDir, 'Fig4'));

%% Fig S1 - behavior (fit and plot psychometric function for each monkey)

% acc = cell(numel(dts), 1);
% dts = io.getDates;
% for ii = 1:numel(dts)
%     stim = io.loadStim(dts{ii}, 'data/stim');
%     ix = stim.goodtrial;
%     xs = stim.dirprob;
%     ys = -(stim.targchosen-2);
%     acc{ii,1} = unique(xs(ix));
%     acc{ii,2} = grpstats(ys(ix), xs(ix));
% end

plot.init;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 2);

weibull = @(x,k,lambda) 1 - exp(-(x/lambda).^k);
lossfcn = @(k,lambda,x,y) sum((y-weibull(x,k,lambda)).^2);

mnkNms = {'Nancy', 'Pat'};
clrs = {[0.2 0.2 0.7], [0.7 0.2 0.2]};
for ii = 1:numel(mnkNms)
    ixm = io.getMonkeyDateFilter(dts, mnkNms(ii));
    Xm = (cell2mat(acc(ixm,1))+1)/2; % map from [-1,1] -> [0,1]
    Ym = cell2mat(acc(ixm,2));
    objfcn = @(th) lossfcn(th(1),th(2),Xm,Ym);
    ths = fmincon(objfcn, 2*[1 1], [], [], [], [], [0 0], [100 100]);
    
    plot(2*Xm-1, Ym, 'o', 'Color', clrs{ii}, 'LineWidth', 2);
    xfine = linspace(min(Xm), max(Xm));
    yhat = weibull(xfine, ths(1), ths(2));    
    plot(2*xfine - 1, yhat, '-', 'Color', clrs{ii}, 'LineWidth', 3);
end
xlabel('Motion Strength');
ylabel('Proportion Right Choices');
set(gca, 'TickDir', 'out');
set(gca, 'XTick', -1:0.5:1);
set(gca, 'YTick', 0:0.25:1);
plot.setPrintSize(gcf, struct('width', 4, 'height', 3.5));
axis square;

%% Fig S2 - load all fits and compare r-sq

% dts = io.getDates('data/fits/space-time');
% fitnms = {'space-time', 'time-only', 'flat'};
% fitstrs = {'ASD', 'ASD', 'Flat'};
% fldsToRm = {'X', 'Xxy'};
% Cells = cell(numel(fitnms),1);
% for ii = 1:numel(fitnms)
%     ccells = tools.makeFitSummaries(['data/fits/' fitnms{ii}], ...
%         [], fitstrs{ii});
%     for jj = 1:numel(fldsToRm)
%         ccells = rmfield(ccells, fldsToRm{jj});
%     end    
%     Cells{ii} = ccells;
% end
% disp('Done');

xi = 2;
yi = 1;

Cx = Cells{xi};
Cy = Cells{yi};
xlbl = fitnms{xi};
ylbl = fitnms{yi};

xs = [Cx.rsq];
ys = [Cy.rsq];

ixCellsToKeep = figure.filterCellsAndPairs(Cy, false);
xs = xs(ixCellsToKeep);
ys = ys(ixCellsToKeep);

plot.init;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 2);
xlim([-0.1 1]);
ylim([-0.1 1]);
plot(xlim, ylim, '-', 'Color', 0.8*ones(3,1), 'LineWidth', 2);
plot(xs, ys, 'k.');
xlabel(['r^2 (' xlbl ')']);
ylabel(['r^2 (' ylbl ')']);
set(gca, 'XTick', 0:0.25:1);
set(gca, 'YTick', 0:0.25:1);
set(gca, 'TickDir', 'out');
plot.setPrintSize(gcf, struct('width', 4, 'height', 3.5));
axis square;

%% Fig S3 - hyperflow

%% Fig S4 - d'/CP as a function of heterogeneity/eccentricity

curCells = cells(ixCellsToKeep);

figS4 = plot.init;

% figS4a = plot.init;
subplot(1,2,1); hold on;
figure.scatterCellFeatures(curCells, 'rf_ecc', 'dPrime');
title('Sensitivity vs. Eccentricity');

% figS4b = figure.scatterCellFeatures(curCells, 'CP', 'rf_ecc');

% figS4b = plot.init;
subplot(1,2,2); hold on;
figure.scatterCellFeatures(curCells, 'rfSpatialVariability', 'dPrime');
title('Sensitivity vs. Heterogeneity');
