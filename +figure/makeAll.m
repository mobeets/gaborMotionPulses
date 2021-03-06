%% init

fitnm = 'space-time-reboot';
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
% fitAllSTRFs('time-only-reboot', false, 'cells ASD time-only', io.getDates, 'MT');
% fitAllSTRFs('flat-reboot', false, 'cells Flat', io.getDates, 'MT');

%% decode with shuffled cell pairs, triplets, quads, etc.

% decoding analysis, with noise corr shuffle
% cellGroups = tools.makeCellGroups(cells, [2 3 4]);

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
% cells = tools.makeFitSummaries(['data/fits/' fitnm]);
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

% add r-squared to pairs
for ii = 1:numel(allPairs)
    cpair = allPairs(ii);
    cell1 = cells(ismember({cells.name}, cpair.cell1));
    cell2 = cells(ismember({cells.name}, cpair.cell2));
    allPairs(ii).minRsq = min([cell1.rsq cell2.rsq]);
end

ixCellsToKeep = figure.filterCellsAndPairs(cells, false, 0.0);
ixPairsToKeep = figure.filterCellsAndPairs(allPairs, true, 0.0);
pairs = allPairs(ixPairsToKeep);

%% summarize experiment, cell, and pair counts

% get all experiments
ixCellsToKeepTmp = figure.filterCellsAndPairs(cells, false);
curCells = cells(ixCellsToKeepTmp);
dts = {curCells.dt};
dts = cellfun(@(x) str2double(x(1:8)), dts);
sessions = unique(dts);
ncells_per_session = grpstats([curCells.index], dts, @numel);
ntrials_per_session = grpstats([curCells.ntrials], dts, @mean);

disp('-----SUMMARY-----');
% summarize experiments, cells, and trials per monkey
isNancy = sessions >= 20150101;
nexps = [sum(~isNancy) sum(isNancy)];
nexps = [nexps sum(nexps)];
fprintf('# sessions: P=%d, N=%d, total=%d\n', nexps);
ncells = [sum(ncells_per_session(~isNancy)) sum(ncells_per_session(isNancy))];
ncells = [ncells sum(ncells)];
fprintf('# cells: P=%d, N=%d, total=%d\n', ncells);

mu = [mean(ncells_per_session(~isNancy)) mean(ncells_per_session(isNancy))];
sd = [std(ncells_per_session(~isNancy)) std(ncells_per_session(isNancy))];
ns = [sqrt(sum(~isNancy)) sqrt(sum(isNancy))];
avgcells = [mu; sd./ns]; avgcells = avgcells(:);
fprintf('avg # cells per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f\n', avgcells);

mu = [mean(ntrials_per_session(~isNancy)) mean(ntrials_per_session(isNancy))];
sd = [std(ntrials_per_session(~isNancy)) std(ntrials_per_session(isNancy))];
ns = [sqrt(sum(~isNancy)) sqrt(sum(isNancy))];
avgtrials = [mu; sd./ns]; avgtrials = avgtrials(:);
fprintf('avg # trials per session: P = %0.1f +/- %0.1f, N = %0.1f +/- %0.1f\n', avgtrials);

% summarize pair filtering as well
ixCellsToKeepInPairs = figure.filterCellsAndPairs(cells, true, 0);
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

doSaveFigs = true;
% cellName = '20140304-MT_12'; % 300-400
cellName = '20150304b-MT_14'; % 600-700; 60-160
startTrial = 600;
endTrial = 700;

% 2b - stimuli
% fig2b = figure.plotStimOrSpikes(cells, cellName, startTrial, endTrial, 'stim');

% 2c - spikes
% fig2c = figure.plotStimOrSpikes(cells, cellName, startTrial, endTrial, 'spikes');

% 2d - all RFs
fig2d = figure.showAllRfs(cells(ixCellsToKeep));

if doSaveFigs
    curSaveDir = fullfile(saveDir, 'Fig2');
    if ~exist(curSaveDir, 'dir')
        mkdir(curSaveDir);
    end
%     export_fig(fig2b, fullfile(curSaveDir, '2b.pdf'));
%     export_fig(fig2c, fullfile(curSaveDir, '2c.pdf'));
    export_fig(fig2d, fullfile(curSaveDir, '2d.pdf'));
end

%% Fig 3 - decoding

% dd = load(decodingFnm);
% cpairs = dd.allPairs;
cpairs = pairs;
doSaveFigs = false;

% 3b
fig3b = figure.deltaDecodingAcc(cpairs);

% 3a, 3c, 3e
xnm = 'rfCorr';
ynm = 'noiseCorrAR';
fig3a = figure.deltaDecodingScatters(cpairs, xnm, ynm, 'both');
fig3c = figure.deltaDecodingScatters(cpairs, xnm, ynm, 'same');
fig3e = figure.deltaDecodingScatters(cpairs, xnm, ynm, 'different');

% 3d, 3f
xnm = 'noiseCorrAR';
ynm = 'scoreGainWithCorrs';
fig3d = figure.deltaDecodingScatters(cpairs, xnm, ynm, 'same');
fig3f = figure.deltaDecodingScatters(cpairs, xnm, ynm, 'different');

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

% beneficial noise corr
A = ([scs.sameTarg] & [scs.noiseCorrAR]<0);
B = (~[scs.sameTarg] & [scs.noiseCorrAR]>0);
disp(['# with beneficial noise corr: ' ppct(sum(A | B))]);
A = ([scs.sameTarg] & [scs.noiseCorrAR]>0);
B = (~[scs.sameTarg] & [scs.noiseCorrAR]<0);
disp(['# with harmful noise corr: ' ppct(sum(A | B))]);

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

% of pairs in diff pools, % with r_RF < 0 and r_sc > 0
scs = pairs;
scs = scs(~[scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] < 0;
B = [scs.noiseCorrAR] > 0;
disp(['# in diff pools with r_{RF} < 0 and r_{sc} > 0: ' ppct(sum(A & B))]);

% of pairs in diff pools, % with r_RF < 0 and r_sc < 0
scs = pairs;
scs = scs(~[scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] < 0;
B = [scs.noiseCorrAR] < 0;
disp(['# in diff pools with r_{RF} < 0 and r_{sc} < 0: ' ppct(sum(A & B))]);

% of pairs in diff pools, % with r_RF > 0 and r_sc > 0
scs = pairs;
scs = scs(~[scs.sameTarg]);
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.rfCorr] > 0;
B = [scs.noiseCorrAR] > 0;
disp(['# in diff pools with r_{RF} > 0 and r_{sc} > 0: ' ppct(sum(A & B))]);

% # with significant \Delta > 0 or \Delta < 0
scs = pairs;
ppct = @(n) sprintf('%d (%0.2f%%) ', [n 100*n/numel(scs)]);
S = [scs.scoreGainWithCorrs];
se = 2*[scs.scoreGainWithCorrs_se];
a = (S - se > 0);
b = (S + se < 0);
disp(['# with significant \Delta > 0: ' ppct(sum(a)) ...
    '; \Delta < 0: ' ppct(sum(b))]);

% correlation between r_sc and r_RF
scs = pairs;
ppct = @(n) sprintf('%d (%0.1f%%) ', [n 100*n/numel(scs)]);
A = [scs.noiseCorrAR];
B = [scs.rfCorr];
ix = ~isnan(A) & ~isnan(B);
[r,pval] = corr(A(ix)', B(ix)');
disp(['correlation between r_{sc} and r_{RF}: ' sprintf('%0.2f', r) ...
    ', p=' sprintf('%0.6f', pval)]);

% avg noise corr
mu = nanmean([pairs.noiseCorrAR]);
sd = nanstd([pairs.noiseCorrAR]);%/sqrt(numel(pairs));
disp(['avg. noise corr: ' sprintf('%0.2f +/- %0.2f', [mu sd])]);

%% Fig 4 - cell pair examples

doSave = false;
figure.plotExamplePairs(pairs, cells, doSave, ...
    fullfile(saveDir, 'Fig4'));

%% Fig S1 - behavior (fit and plot psychometric function for each monkey)

% dts = io.getDates;
% acc = cell(numel(dts), 4);
% for ii = 1:numel(dts)
%     stim = io.loadStim(dts{ii}, 'data/stim');
%     ix = stim.goodtrial;
%     xs = stim.dirprob;
%     
%     Xs = sum(stim.pulses(ix,:,:),3);
%     sd = std(Xs(:)); % standard deviation over all pulses
%     X = Xs/sd;
%     coh = mean(X,2);
%     
%     ys = -(stim.targchosen-2);
%     acc{ii,1} = unique(xs(ix));
%     acc{ii,2} = grpstats(ys(ix), xs(ix));
%     acc{ii,3} = coh;
%     acc{ii,4} = ys(ix);
% end

% From Yates2017:
% The strength of the pulses was normalized by the s.d. of all 
% pulse values shown. The psychophysical performance (Fig. 1d) was
% measured by calculating the proportion of rightward choices as a
% function of the net motion on each trial (sum of normalized pulse
% strength). The net motion strengths from all trials were divided
% into 30 equal quantiles to form the bins in Figure 1d.

saveFig = false;
plot.init;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 2);

weibull = @(x,k,lambda) 1 - exp(-(x/lambda).^k);
lossfcn = @(k,lambda,x,y) sum((y-weibull(x,k,lambda)).^2);
doFit = false;

if doFit
    xind = 1; yind = 2;
else
    xind = 3; yind = 4;
end
mnkNms = {'Nancy', 'Pat', ''};
clrs = {[0.2 0.2 0.7], [0.7 0.2 0.2], 0.2*ones(1,3)};
for ii = 3%1:2%numel(mnkNms)
    if ~isempty(mnkNms{ii})
        ixm = io.getMonkeyDateFilter(dts, mnkNms(ii));
    else
        ixm = true(size(acc,1),1);
    end
    Xm = cell2mat(acc(ixm,xind));
    Ym = cell2mat(acc(ixm,yind));
    
    if doFit
        Xm = (Xm+1)/2; % map from [-1,1] -> [0,1]
        objfcn = @(th) lossfcn(th(1),th(2),Xm,Ym);
        ths = fmincon(objfcn, 1*[1 1], [], [], [], [], [0 0], [100 100]);
        Xc = 2*Xm-1;
        Yc = Ym;
    else
        nbins = 30;
        bins = prctile(Xm, linspace(0, 100, nbins));
%         bins = linspace(min(Xm), max(Xm), nbins);% bins(end) = bins(end)+0.1;
%         bins = linspace(prctile(Xm,0.5), prctile(Xm,99.5), nbins);
        bins = [bins inf];
        Yc = nan(nbins,1);
        ses = nan(nbins,1);
        ns = nan(nbins,1);
        for jj = 1:(numel(bins)-1)
            ixc = (Xm >= bins(jj)) & (Xm < bins(jj+1));
            Yc(jj) = nanmean(Ym(ixc));
            ses(jj) = nanstd(Ym(ixc))/sqrt(sum(ixc));
            ns(jj) = sum(ixc);
        end
        Xc = bins(1:end-1);
    end    
    plot(Xc, Yc, 'o', 'MarkerFaceColor', clrs{ii}, 'Color', clrs{ii}, ...
        'LineWidth', 2);
    if ~doFit
        for jj = 1:numel(ses)
            plot(Xc(jj)*[1 1], Yc(jj) + ses(jj)*[-1 1], 'k-', ...
                'LineWidth', 2, 'Color', clrs{ii});
        end

        [logitCoef,dev] = glmfit(Xm, Ym, 'binomial', 'logit');
        logitFit = glmval(logitCoef, Xc', 'logit');
        plot(Xc, logitFit, '-', 'Color', clrs{ii}, 'LineWidth', 3);
    end
    
    if doFit
        xfine = linspace(min(Xm), max(Xm));
        yhat = weibull(xfine, ths(1), ths(2));
        Xcfine = 2*xfine - 1;
        plot(Xcfine, yhat, '-', 'Color', clrs{ii}, 'LineWidth', 3);
    end
end
xlabel('cumulative motion strength');
ylabel('proportion right choices');
set(gca, 'TickDir', 'out');
if doFit
    set(gca, 'XTick', -1:0.5:1);
else
    set(gca, 'XTick', -2:1:2);
    axis tight;
end
set(gca, 'YTick', 0:0.25:1);
plot.setPrintSize(gcf, struct('width', 4, 'height', 3.5));
axis square;

if saveFig
    fnm = fullfile(saveDir, 'Fig1D.pdf');
    export_fig(gcf, fnm);
end

%% Fig S2 - load all fits and compare r-sq

% dts = io.getDates('data/fits/space-time-reboot');
% % fitnms = {'space-time', 'time-only', 'flat'};
% fitnms = {'space-time-reboot', 'time-only-reboot', 'flat-reboot'};
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

doSaveFigs = false;
xi = 2;
yi = 1;

Cx = Cells{xi};
Cy = Cells{yi};

lblnms = {'spatial+temporal RF', 'temporal RF', 'flat'};
xlbl = lblnms{xi};
ylbl = lblnms{yi};

xs = [Cx.rsq];
ys = [Cy.rsq];

ixCellsToKeep = figure.filterCellsAndPairs(Cy, false);
xs = xs(ixCellsToKeep);
ys = ys(ixCellsToKeep);

figSx = plot.init;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 2);
xlim([-0.1 1]); ylim(xlim);
plot(xlim, ylim, '-', 'Color', 0.8*ones(3,1), 'LineWidth', 2);
plot(xs, ys, 'k.', 'MarkerSize', 10);
xlabel(['r^2 (' xlbl ')']);
ylabel(['r^2 (' ylbl ')']);
set(gca, 'XTick', 0:0.25:1);
set(gca, 'YTick', 0:0.25:1);
set(gca, 'TickDir', 'out');
plot.setPrintSize(gcf, struct('width', 4, 'height', 3.5));
axis square;
if doSaveFigs
    fnm = fullfile(saveDir, 'FigSx_r2.pdf');
    export_fig(figSx, fnm);
end

%% Fig S3 - hyperflow with arrows

doSaveFigs = false;

exampleCellNames = {'20140304-MT_6', '20150324a-MT_14', '20150324a-MT_13', ...
    '20140304-MT_2', '20140307-MT_1', '20140304-MT_3', '20150324a-MT_5'};
% exampleCellNames = {'20150306b-MT_9'};
exampleCellNames = {'20140303-MT_3', '20140305-MT_2', '20140305-MT_4', ...
    '20140306-MT_2', '20150324a-MT_4', '20150324a-MT_10', '20150324a-MT_12'};
exampleCellNames = {'20140213-MT_4', '20140310-MT_5'};

for kCell = 1:numel(exampleCellNames)
    exampleCellMT = exampleCellNames{kCell};
    vMT = cells(strcmp({cells.name}, exampleCellMT));
    if isempty(vMT)
        continue;
    end
    d = io.loadDataByDate(vMT.dt);
    n = d.neurons{vMT.index};
    if isempty(n.hyperflow)
        continue
    end
    t1 = d.stim.targ1XY; t2 = d.stim.targ2XY;
    
    figure;
%     set(gcf, 'color', 'w');
    plot.getColors([0 1]);
    plot.plotHyperflowMT(n, t1, t2, 5); axis xy;
    xd = xlim; yd = ylim;
    axis equal; xlim(xd); ylim(yd); axis off;
    figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [50 50]);
    if doSaveFigs
        fignm = vMT.name;
        cSaveDir = fullfile(saveDir, 'hyperflow2');
        if ~exist(cSaveDir, 'dir')
            mkdir(cSaveDir);
        end
        fnm = fullfile(cSaveDir, [fignm '.pdf']);
        export_fig(fnm);
        close;
    end
end

% todo:
% - need to make sure the pdfs in the .ai file are the ones in that folder
% - find cells with as much subfields as possible
% - confirm sign of hyperflow vs. RF (e.g., 20140310-MT_5)

%% find all cells with hyperflow

% ns = dir('data/neurons/');
% Ns = [];
% for ii = 5:numel(ns)
%     n = load(['data/neurons/' ns(ii).name]);
%     if isfield(n, 'hyperflow') && ~isempty(n.hyperflow)
%         Ns = [Ns n];
%     end
% end
% dtsWithHyperflow = unique({Ns.exname});
exampleCellNames = {};
% goodCells = cells([cells.rsq] < 0.0);
goodCells = cells(ixCellsToKeep);
for ii = 1:numel(dtsWithHyperflow)
    cCells = {goodCells(ismember({goodCells.dt}, dtsWithHyperflow{ii}(2:end))).name};
    exampleCellNames = [exampleCellNames cCells];
end

%% Fig S3 - hyperflow with RF

doSaveFigs = false;
% exampleCellNames = {'20150324a-MT_14'};
% exampleCellNames = {'20150518-MT_12'};
% exampleCellNames = {'20140310-MT_5'};

% ids = [Ns(ismember({Ns.exname}, cdt)).id];

showTargs = false;
showHyperflow = true;
for kCell = 1:numel(exampleCellNames)
    exampleCellMT = exampleCellNames{kCell};
    vMT = cells(strcmp({cells.name}, exampleCellMT));
    if isempty(vMT)
        continue;
    end
    vMT.label = 'MT'; vMT.score = vMT.rsq;
    
    % hyperflow ASD overlay
    d = io.loadDataByDate(vMT.dt);
    n = d.neurons{vMT.index};
    figure;    
    plot.plotSaccadeKernelOverlay(d.stim, n, vMT, showTargs, ...
        showHyperflow, 1);
    title('');
    xd = xlim; yd = ylim;
    xlim(xd); ylim(yd); axis off;
    figure.cleanupForPrint(gcf, 'FontSize', 8, 'PaperSize', [50 50])
    if doSaveFigs
        fignm = vMT.name;
        cSaveDir = fullfile(saveDir, 'hyperflow2');
        if ~exist(cSaveDir, 'dir')
            mkdir(cSaveDir);
        end
        fnm = fullfile(cSaveDir, [fignm '_withRF.pdf']);
        export_fig(fnm);
        close;
     end
end

%% Fig S4 - d'/CP as a function of heterogeneity/eccentricity

doSaveFigs = true;
fnma = fullfile(saveDir, 'FigS3_a.pdf');
fnmb = fullfile(saveDir, 'FigS3_b.pdf');

curCells = cells(ixCellsToKeep);
% figS4 = plot.init;

figS4a = plot.init;
% subplot(1,2,1); hold on;
figure.scatterCellFeatures(curCells, 'rf_ecc', 'dPrime');
axis tight; yl = ylim; ylim([0 yl(2)]);
% title('Sensitivity vs. Eccentricity');
if doSaveFigs
    export_fig(figS4a, fnma);
end

% figS4b = figure.scatterCellFeatures(curCells, 'CP', 'rf_ecc');

figS4b = plot.init;
% subplot(1,2,2); hold on;
figure.scatterCellFeatures(curCells, 'rfSpatialVariability', 'dPrime');
axis tight; yl = ylim; ylim([0 yl(2)]);
% title('Sensitivity vs. Heterogeneity');
if doSaveFigs
    export_fig(figS4b, fnmb);
end
