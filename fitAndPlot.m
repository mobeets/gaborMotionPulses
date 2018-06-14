%% load cell fits

fitnm = 'new-2';
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');
cells = tools.makeFitSummaries(['data/fits/' fitnm]);

%% count cells

dts = {cells.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);

ix = ([cells.pctCorrect] >= 0.70);
ix = ix & ([cells.ntrials] >= 100);
ix = ix & (abs([cells.dPrime]) >= 0.1);

sessions = unique(dts(ix));
ncells = grpstats([cells(ix).index], dts(ix), @numel);
ntrials = grpstats([cells(ix).ntrials], dts(ix), @mean);
[sessions' ncells]

ixIsNancy = sessions >= 20150101;
[sum(~ixIsNancy) sum(ixIsNancy)]
[sum(ncells(~ixIsNancy)) sum(ncells(ixIsNancy))]
[mean(ncells(~ixIsNancy)) mean(ncells(ixIsNancy))]
[std(ncells(~ixIsNancy)) std(ncells(ixIsNancy))]
[mean(ntrials(~ixIsNancy)) mean(ntrials(ixIsNancy))]

isSpatialVar = log([cells.rfSpatialVariability]) > -15;
ix = ix & isSpatialVar;
[sum(ix & dts < 20150101) sum(ix & dts >= 20150101)]

%% find cell pairs

allPairs2 = tools.makeCellPairs(cells);
for ii = 1:numel(allPairs)
    allPairs(ii).rfCorr_spatial = allPairs2(ii).rfCorr_spatial;
    allPairs(ii).sameCorr_spatial = allPairs2(ii).sameCorr_spatial;
end

%% compute delta decoding accuracy

nShuffles = 25;
scs = tools.decodeAndShuffle({allPairs.stimdir}, ...
    {allPairs.Ys}, nShuffles);
% [allScsRaw, allScsShuf] = tools.decodeAndShuffle({allPairs.stimdir}, ...
%    {allPairs.Ys}, nShuffles, {allPairs.dirstrength_binned});

% save('data/decode/decodeAndShuffle.mat', 'scs', 'allPairs');
load('data/decode/decodeAndShuffle.mat', 'scs', 'allPairs');

scDelta = num2cell(scs.scsDelta);
[allPairs.scoreGainWithCorrs] = scDelta{:};
scDeltaSe = num2cell(scs.scsShufSe);
[allPairs.scoreGainWithCorrs_se] = scDeltaSe{:};

%% filtering pairs

ix = fig.filterPairs(allPairs);
pairs = allPairs(ix);

dts = {pairs.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
[sum(dts < 20150101) sum(dts >= 20150101)]

%% plot score change with shuffles

doSave = false;
saveDir = 'data/figs/summary/decoding';
if doSave && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

fig.deltaDecodingAcc; % must do this before the below

fnm = fullfile(saveDir, 'scores.pdf');
export_fig(gcf, fnm);

%% plot scatters

doSave = false;
saveDir = 'data/figs/summary/decoding';
if doSave && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

xnms = {'rfCorr', 'noiseCorrAR'};
% xnms = {'rfCorr_spatial', 'noiseCorrAR'};
ynms = {'noiseCorrAR', 'scoreGainWithCorrs'};
fnms = {'rscVsRf', 'scoreVsrsc'};
for ii = 1%1:numel(xnms)
    xnm = xnms{ii};
    ynm = ynms{ii};
    fnm = fnms{ii};
    fig.deltaDecodingScatters;
    
    if doSave
        fnm1 = fullfile(saveDir, [fnm '-diffPool.pdf']);
        fnm2 = fullfile(saveDir, [fnm '-samePool.pdf']);
        fnm3 = fullfile(saveDir, [fnm '-all.pdf']);
        export_fig(fig1, fnm1);
        export_fig(fig2, fnm2);
        export_fig(fig3, fnm3);
    end
end

%% d' and CP vs. RF center

fig.DprimeAndCp;

%% plot example cells where correlations helped

fig.plotExamplePairs;
