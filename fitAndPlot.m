%% load cell fits

fitnm = 'new-2';
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');
cells = tools.makeFitSummaries(['data/fits/' fitnm]);

%% find cell pairs

allPairs = tools.makeCellPairs(cells);

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

doSave = true;
saveDir = 'data/figs/summary/decoding';
if doSave && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

xnms = {'rfCorr', 'noiseCorrAR'};
ynms = {'noiseCorrAR', 'scoreGainWithCorrs'};
fnms = {'rscVsRf', 'scoreVsrsc'};
for ii = 2%1:numel(xnms)
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
