%% run fits

fitnm = 'new-2';
% fitAllSTRFs(fitnm, false, 'cells ASD', io.getDates, 'MT');

%% load fits

fitdir = fullfile('data', 'fits', fitnm);
vu = tools.makeFitSummaries(fitdir, true, 'ASD');
vu = figure.filterData(vu);

%% only keep cells with spatially varying RFs

tm = arrayfun(@(v) var(v.wfSvd_1(:,1)), vu, 'uni', 0);
[vu.wfVar] = tm{:};
vuMT = vu([vu.isMT]);
vuMT_var = vuMT(log([vuMT.wfVar])>-17);

%% decoding

nShuffles = 20;
pairs = tools.makeCellPairs(vuMT_var);
[scsRaw, scsShuf] = tools.decodeAndShuffle({pairs.stim}, {pairs.Ys}, nShuffles);

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
plot(scsDelta(ix,:));
plot(xlim, [0 0], 'k-');
xlabel('pair index (sorted)');
ylabel('decoding improvement compared to shuffled');

%% plot

doSave = false;
scs = pairs;

ix = [scs.sameTarg];
ix1 = [scs.sameCorr];
ix2a = (ix&~ix1); % same dir, anti-corr
ix2b = (~ix&ix1); % opp dir, pos-cor
ix2 = ix2a | ix2b;

xnms = {'rfCorr', 'noiseCorrAR', 'rfCorr', 'rfDist'};
ynms = {'noiseCorrAR', 'scoreGainWithCorrs', 'scoreGainWithCorrs', 'scoreGainWithCorrs'};
xlbls = {'RF correlation', 'noise correlation', 'RF correlation', 'RF center dist'};
ylbls = {'noise correlation', '\Delta decoding accuracy', '\Delta decoding accuracy', '\Delta decoding accuracy'};
fnms = {'rscVsRf', 'scoreVsrsc', 'scoreVsRf', 'rfdistVsScore'};

for ii = 1:numel(xnms)
    xnm = xnms{ii};
    ynm = ynms{ii};
    xlbl = xlbls{ii};
    ylbl = ylbls{ii};
    fnm = fnms{ii};
    fnm = fullfile('data', 'figs', 'summary', [fitnm '_' fnm '.pdf']);

    figure; hold on; set(gcf, 'color', 'w');
    t1 = scs(ix2a); % same dir, anti-corr
    t2 = scs(ix2b); % opp dir, pos-cor
    t3 = scs(~ix2 & ~ix);
    t4 = scs(~ix2 & ix);
    c2 = [0.3 0.6 0.6];
    c1 = [0.9 0.7 0.2];

    nms = {'different pools, -RF corr', 'different pools, +RF corr', ...
        'same pool, -RF corr', 'same pool, +RF corr'};

    xmx = max(abs([scs.(xnm)]));
    ymx = max(abs([scs.(ynm)]));
    xl = [-xmx xmx];
    yl = [-ymx ymx];

    subplot(2,2,1); hold on; set(gca, 'FontSize', 14);
    scatter([t3.(xnm)], [t3.(ynm)], 50, c1, 'filled');
    xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
    xlim(xl); ylim(yl);
    title(nms{1});

    subplot(2,2,2); hold on; set(gca, 'FontSize', 14);
    scatter([t2.(xnm)], [t2.(ynm)], 50, c1);
    xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
    xlim(xl); ylim(yl);
    title(nms{2});

    subplot(2,2,3); hold on; set(gca, 'FontSize', 14);
    scatter([t1.(xnm)], [t1.(ynm)], 50, c2);
    xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
    xlim(xl); ylim(yl);
    title(nms{3});

    subplot(2,2,4); hold on; set(gca, 'FontSize', 14);
    scatter([t4.(xnm)], [t4.(ynm)], 50, c2, 'filled');
    xlabel(xlbl); ylabel(ylbl); plot(xl, [0 0], 'k--'); plot([0 0], yl, 'k--');
    xlim(xl); ylim(yl);
    title(nms{4});
    
    if doSave
        export_fig(gcf, fnm);
    end
end
