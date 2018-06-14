%% init

fitnm = 'new-2';
doSaveFigs = true;

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
disp(sprintf('# sessions: P=%d, N=%d', nexps));
ncells = [sum(ncells_per_session(~isNancy)) sum(ncells_per_session(isNancy))];
disp(sprintf('# cells: P=%d, N=%d', ncells));

mu = [mean(ncells_per_session(~isNancy)) mean(ncells_per_session(isNancy))];
sd = [std(ncells_per_session(~isNancy)) std(ncells_per_session(isNancy))];
avgcells = [mu; sd]; avgcells = avgcells(:);
disp(sprintf('avg # cells per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f', avgcells));

mu = [mean(ntrials_per_session(~isNancy)) mean(ntrials_per_session(isNancy))];
sd = [std(ntrials_per_session(~isNancy)) std(ntrials_per_session(isNancy))];
avgtrials = [mu; sd]; avgtrials = avgtrials(:);
disp(sprintf('avg # trials per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f', avgtrials));

% summarize pair filtering as well
ixCellsToKeepInPairs = figure.filterCellsAndPairs(cells, true);
curCells = cells(ixCellsToKeepInPairs);
dts = {curCells.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
isNancy = dts >= 20150101;
ncells = [sum(~isNancy) sum(isNancy)];
disp(sprintf('# cells used in pairs: P=%d, N=%d', ncells));

dts = {pairs.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
isNancy = dts >= 20150101;
ncells = [sum(~isNancy) sum(isNancy)];
disp(sprintf('# pairs: P=%d, N=%d', ncells));

%% Fig 2 - ASD

cellName = '20140304-MT_12';
startTrial = 300;
endTrial = 400;

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

%% Fig 4 - cell pair examples

figure.plotExamplePairs(pairs, cells, doSaveFigs, ...
    fullfile(saveDir, 'Fig4'));

%% Fig Sx - hyperflow

%% Fig Sx - CP, dPrime, etc.

figS1 = figure.dPrimeAndCP(cells, 'dPrime');
figS2 = figure.dPrimeAndCP(cells, 'CP');
