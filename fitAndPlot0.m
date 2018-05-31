%% load

fitnm = 'new-2';
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');
cells = tools.makeFitSummaries(fitdir);
% ix = ix & ([cells.pctCorrect] >= 0.75);
% ix = ix & ([cells.ntrials] >= 100);
% ix = ix & ([cells.dPrime] >= 0);

%% run fits

fitnm = 'new-2';
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');

%% load fits

fitdir = fullfile('data', 'fits', fitnm);
vu = tools.makeFitSummaries(fitdir, true, 'ASD');
vu = figure.filterData(vu);

%% only keep cells with spatially varying RFs

% todo: should use wfSvd_U(:,1)
% wfSvd_1 is just the rank-1 RF for the first temporal pulse...
% tm = arrayfun(@(v) var(v.wfSvd_1(:,1)), vu, 'uni', 0);
tm = arrayfun(@(v) var(v.wfSvd_U(:,1)), vu, 'uni', 0);
[vu.wfVar] = tm{:};
vuMT = vu([vu.isMT]);
vuMT_var = vuMT(log([vuMT.wfVar])>-17);

%% decoding

nShuffles = 10;
pairs = tools.makeCellPairs(vuMT_var);
[scsRaw, scsShuf] = tools.decodeAndShuffle({pairs.stimdir}, {pairs.Ys}, ...
    nShuffles);
% [scsRaw, scsShuf] = tools.decodeAndShuffle({pairs.stimdir}, {pairs.Ys}, ...
%     nShuffles, {pairs.dirstrength_binned});

%% check noise corr

rscs = nan(numel(pairs), 3);
for ii = 1:numel(pairs)
    x = pairs(ii).stimdir;
    ys = pairs(ii).Ys;
    ix = ~any(isnan(x),2);
    x = x(ix,:); ys = ys(ix,:);
    
    r1 = corr(ys(x == -1,:)); r1 = r1(2);
    r2 = corr(ys(x == 1,:)); r2 = r2(2);
    r3 = pairs(ii).noiseCorrAR;
    rscs(ii,:) = [r1 r2 r3];
    pairs(ii).rscL = r1;
    pairs(ii).rscR = r2;
    pairs(ii).rscMatch = (sign(r1) == sign(r3)) & (sign(r2) == sign(r3));
end

%% plot score change with shuffles

% compute score change from shuffle
scsShufMean = mean(scsShuf, 2);
scsShufSe = (std(scsShuf,[],2)/sqrt(size(scsShuf,2)));
scsDelta = bsxfun(@plus, -[(scsShufMean+2*scsShufSe) scsShufMean ...
    (scsShufMean-2*scsShufSe)], scsRaw);
scsDeltaMean = num2cell(scsDelta(:,2));
[pairs.scoreGainWithCorrs] = scsDeltaMean{:};

figure; set(gcf, 'color', 'w'); set(gca, 'FontSize', 16); hold on;
[~,ix] = sort(scsDelta(:,2));
plot([1 size(scsDelta,1)], [0 0], '-', 'Color', 0.8*ones(3,1));
plot(scsDelta(ix,:), '-', 'Color', [0.4 0.4 0.4]);
xlabel('pair index (sorted)');
ylabel({'\Delta decoding accuracy', ...
    '\leftarrow Correlations hurt        Correlations help \rightarrow'});

%% plot scatters

scs = pairs;
doSave = false;

% ix = ([scs.cell1_dPrime] > 0.5) & ([scs.cell2_dPrime] > 0.5);
% ix = [scs.rscMatch];
% scs = scs(ix);

samePool = [scs.sameTarg];
posRfCorr = [scs.sameCorr];
ix1 = ~samePool & ~posRfCorr;
ix2 = ~samePool & posRfCorr;
ix3 = samePool & ~posRfCorr;
ix4 = samePool & posRfCorr;

xnms = {'rfCorr', 'noiseCorrAR', 'rfCorr', 'rfDist'};
ynms = {'noiseCorrAR', 'scoreGainWithCorrs', 'scoreGainWithCorrs', 'noiseCorrAR'};
fnms = {'rscVsRf', 'scoreVsrsc', 'scoreVsRf', 'rfdistVsrsc'};
cnm = 'rfCorr';

% xnms = {'rfDist'};
% ynms = {'noiseCorrAR'};
% fnms = {'rscVsRf'};

fldnms = struct('rfCorr', 'RF correlation', ...
    'noiseCorrAR', 'noise correlation', ...
    'rfDist', 'RF center distance', ...
    'scoreGainWithCorrs', '\Delta decoding accuracy');
fldnms = containers.Map(fieldnames(fldnms), struct2cell(fldnms));

mkrsz = 50;
nclrs = 21;
cs1 = cbrewer('seq', 'Reds', nclrs);
cs2 = cbrewer('seq', 'Blues', nclrs);
valToClrInd = @(v, vmn, vmx) floor((nclrs-1)*(v - vmn)/(vmx - vmn))+1;

for ii = 1:numel(xnms)
    xnm = xnms{ii};
    ynm = ynms{ii};
    xlbl = fldnms(xnm);
    ylbl = fldnms(ynm);
    fnm = fnms{ii};
    fnm = fullfile('data', 'figs', 'summary', [fitnm '_' fnm '.png']);

    figure; hold on; set(gcf, 'color', 'w');
    t1 = scs(ix1);
    t2 = scs(ix2);
    t3 = scs(ix3);
    t4 = scs(ix4);
%     c1 = [0.8 0.3 0.3];
%     c2 = [0.3 0.3 0.8];    

    xmx = max(abs([scs.(xnm)]));
    ymx = max(abs([scs.(ynm)]));
    xl = [-xmx xmx];
    yl = [-ymx ymx];
    
    nms = {'different pools', 'same pool'};
    nrows = 2; ncols = 1;

    subplot(nrows,ncols,1); hold on; set(gca, 'FontSize', 14);
    clrs = cs1(valToClrInd([t1.(cnm)], -1, 1),:);
    scatter([t1.(xnm)], [t1.(ynm)], mkrsz, clrs, 'filled');
    clrs = cs1(valToClrInd([t2.(cnm)], -1, 1),:);
    scatter([t2.(xnm)], [t2.(ynm)], mkrsz, clrs, 'filled');
    xlabel(xlbl); ylabel(ylbl);
    plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
    plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
    xlim(xl); ylim(yl);
    set(gca, 'TickLength', [0 0]);
    title(nms{1});

    subplot(nrows,ncols,2); hold on; set(gca, 'FontSize', 14);
    clrs = cs2(valToClrInd([t3.(cnm)], -1, 1),:);
    scatter([t3.(xnm)], [t3.(ynm)], mkrsz, clrs, 'filled');
    clrs = cs2(valToClrInd([t4.(cnm)], -1, 1),:);
    scatter([t4.(xnm)], [t4.(ynm)], mkrsz, clrs, 'filled');
    xlabel(xlbl); ylabel(ylbl);
    plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
    plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
    xlim(xl); ylim(yl);
    set(gca, 'TickLength', [0 0]);
    title(nms{2});    
    
    plot.setPrintSize(gcf, struct('width', 4, 'height', 6));
    
%     nms = {'different pools, -RF corr', 'different pools, +RF corr', ...
%         'same pool, -RF corr', 'same pool, +RF corr'};    
%     nrows = 2; ncols = 2;

%     subplot(nrows,ncols,1); hold on; set(gca, 'FontSize', 14);
%     scatter([t1.(xnm)], [t1.(ynm)], 50, c1, 'filled');
%     xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
%     xlim(xl); ylim(yl);
%     title(nms{1});
% 
%     subplot(nrows,ncols,2); hold on; set(gca, 'FontSize', 14);
%     scatter([t2.(xnm)], [t2.(ynm)], 50, c1);
%     xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
%     xlim(xl); ylim(yl);
%     title(nms{2});
% 
%     subplot(nrows,ncols,3); hold on; set(gca, 'FontSize', 14);
%     scatter([t3.(xnm)], [t3.(ynm)], 50, c2);
%     xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
%     xlim(xl); ylim(yl);
%     title(nms{3});
% 
%     subplot(nrows,ncols,4); hold on; set(gca, 'FontSize', 14);
%     scatter([t4.(xnm)], [t4.(ynm)], 50, c2, 'filled');
%     xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
%     xlim(xl); ylim(yl);
%     title(nms{4});
    
    if doSave
        export_fig(gcf, fnm);
    end
end

%% d' and CP vs. RF center

doSave = false;
fnm = fullfile('data', 'figs', 'summary', [fitnm '_dpAndCp.pdf']);

nrows = 1; ncols = 3;
figure; hold on; set(gcf, 'color', 'w');

subplot(nrows,ncols,1); hold on; set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);
xs = [vuMT_var.rf_ecc];
ys = [vuMT_var.dPrime];
mdl = fitlm(xs, ys);
h = mdl.plot;
h(1).Marker = 'o';
h(1).Color = 'k';
h(1).MarkerSize = 4.5;
h(1).MarkerFaceColor = 0.5*ones(3,1);
h(1).LineWidth = 1.5;
h(2).LineWidth = 2;
h(2).Color = [0.8 0.2 0.2];
set(gca, 'TickDir', 'out');
legend off;
xlabel('RF eccentricity');
ylabel('d''');
xlim([0 1.1]);
title('');

subplot(nrows,ncols,2); hold on; set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);
xs = [vuMT_var.rf_ecc];
ys = [vuMT_var.cp_YresARc];
xs = xs(ys > 0); ys = ys(ys > 0);
plot([min(xs) max(xs)], [0.5 0.5], 'k--', 'LineWidth', 2);
mdl = fitlm(xs, ys);
h = mdl.plot;
h(1).Marker = 'o';
h(1).Color = 'k';
h(1).MarkerSize = 4.5;
h(1).MarkerFaceColor = 0.5*ones(3,1);
h(1).LineWidth = 1.5;
h(2).LineWidth = 3;
h(2).Color = [0.8 0.2 0.2];
set(gca, 'TickDir', 'out');
legend off;
xlabel('RF eccentricity');
ylabel('CP');
xlim([0 1.1]);
% ylim([0.4 0.7]);
title('');

subplot(nrows,ncols,3); hold on; set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);
xs = [pairs.rfCorr];
ys = [pairs.noiseCorrAR];
plot([min(xs) max(xs)], [0 0], 'k--', 'LineWidth', 2);
mdl = fitlm(xs, ys);
h = mdl.plot;
h(1).Marker = 'o';
h(1).Color = 'k';
h(1).MarkerSize = 4.5;
h(1).MarkerFaceColor = 0.5*ones(3,1);
h(1).LineWidth = 1.5;
h(2).LineWidth = 3;
h(2).Color = [0.8 0.2 0.2];
set(gca, 'TickDir', 'out');
legend off;
xlabel('RF correlation');
ylabel('Noise correlation');
title('');

plot.setPrintSize(gcf, struct('width', 11, 'height', 3));

if doSave
    export_fig(gcf, fnm);
end

%% plot example cells where correlations helped

doSave = false;
fnm = fullfile('data', 'figs', 'summary', [fitnm '_examples.pdf']);

scs = pairs;
ix = ([scs.cell1_dPrime] > 0.5) & ([scs.cell2_dPrime] > 0.5);
% ix = [scs.rscMatch];
scs = scs(ix);

samePool = [scs.sameTarg];
posNoiseCorr = [scs.noiseCorrAR] > 0;
posDecAcc = [scs.scoreGainWithCorrs] > 0.02;

% diff pools, beneficial noise corrs
ixA = ~samePool & posNoiseCorr & posDecAcc;
psA = scs(ixA);
psA = psA([2 4]);

% diff pools, harmful noise corrs
ixB = ~samePool & ~posNoiseCorr & [scs.scoreGainWithCorrs] < -0.02;
psB = scs(ixB);

% same pools, harmful noise corrs
ixC = samePool & posNoiseCorr & [scs.scoreGainWithCorrs] < -0.05;
psC = scs(ixC);
psC = psC(1);

% same pools, helpful noise corrs
ixD = samePool & posNoiseCorr & [scs.scoreGainWithCorrs] > 0.02;
psD = scs(ixD);

ps = [psC psB psA];
% ps = psD;
nms = {'Same pool, +r_{sc}', 'Different pools, -r_{sc}', ...
    'Different pools, +r_{sc}', 'Different pools, +r_{sc}'};

figure; set(gcf, 'color', 'w');
for ii = 1:numel(ps)
    subplot(2, ceil(numel(ps)/2), ii); hold on; set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 2);
    S = plot.visualizePairwiseCorr(ps(ii).Ys, ps(ii).stimdir, false);
    sc1 = ps(ii).scoreGainWithCorrs;
    sc2 = ps(ii).noiseCorrAR;
    sc3 = ps(ii).rfCorr;
    nm = [num2str(ps(ii).cell1) ', ' num2str(ps(ii).cell2)];
    vals = [sprintf('r_{sc} = %0.2f', sc2) ...
        ', ' sprintf('r_{RF} = %0.2f', sc3) ...
        ', \Delta acc = ' sprintf('%0.1f%%', 100*sc1)];
%     title({nm, vals});
    title(vals);
    set(gca, 'TickDir', 'out');
    xlabel('Neuron 1 spike count');
    ylabel('Neuron 2 spike count');
end
plot.setPrintSize(gcf, struct('width', 7, 'height', 6));

if doSave
    export_fig(gcf, fnm);
end

%%

ps = psA;

for ii = 1:numel(ps)
    figure;
    subplot(2,2,4); hold on;
    ix1 = strcmp({vuMT_var.name}, ps(ii).cell1);
    v = vuMT_var(ix1);
    plot.plotCellRF_fit(v.Xxy, v.wfSvd_1(:,1), nan, 20);
    subplot(2,2,1); hold on;
    ix1 = strcmp({vuMT_var.name}, ps(ii).cell2);
    v = vuMT_var(ix1);
    plot.plotCellRF_fit(v.Xxy, v.wfSvd_1(:,1), nan, 20);
    
    subplot(2,2,2); hold on; set(gca, 'FontSize', 16);
    set(gca, 'LineWidth', 2);
    S = plot.visualizePairwiseCorr(ps(ii).Ys, ps(ii).stimdir, false);
    sc1 = ps(ii).scoreGainWithCorrs;
    sc2 = ps(ii).noiseCorrAR;
    sc3 = ps(ii).rfCorr;
    nm = [num2str(ps(ii).cell1) ', ' num2str(ps(ii).cell2)];
    vals = [sprintf('r_{sc} = %0.2f', sc2) ...
        ', ' sprintf('r_{RF} = %0.2f', sc3) ...
        ', \Delta acc = ' sprintf('%0.1f%%', 100*sc1)];
%     title({nm, vals});
    title(vals);
    set(gca, 'TickDir', 'out');
    xlabel('Neuron 1 spike count');
    ylabel('Neuron 2 spike count');
    plot.setPrintSize(gcf, struct('width', 7, 'height', 6));
end
