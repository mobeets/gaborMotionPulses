%% init

fitnm = 'new-2';
doSaveFigs = false;

saveDir = fullfile('data', 'figures', fitnm);
if doSaveFigs && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
pairsFnm = fullfile(saveDir, 'allPairs.mat');
decodingFnm = fullfile(saveDir, 'allPairsWithScs.mat');
% addpath(genpath('~/code/gaborMotionASD/mASD/'));

%% make fits, cell pairs, and perform decoding analysis

% fits
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');

% cell pairs
% allPairs = tools.makeCellPairs(cells);
% save(pairsFnm, 'allPairs');

% decoding analysis, with noise corr shuffle
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


%% Fig 4 - cell pair examples

figure.plotExamplePairs(pairs, cells, doSaveFigs, ...
    fullfile(saveDir, 'Fig4'));

%% Fig Sx - hyperflow

%% Fig Sx - CP, dPrime, etc.

figS1 = figure.dPrimeAndCP(cells, 'dPrime');
figS2 = figure.dPrimeAndCP(cells, 'CP');
