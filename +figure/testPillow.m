%% load cell pairs

% pairs_all = tools.makeCellPairs(cells);
pairs_all = allPairs;
ixPairsToKeep = figure.filterCellsAndPairs(pairs_all, true, 0);
pairs = pairs_all(ixPairsToKeep);
dts = unique({pairs.dt});

%% part 1. RFs predict noise correlation

% current status: I think this is done. we don't need to do triplets or
% more (i.e., pairs is fine) because the claim is only that the dot-product
% of the RF maps is proportional to the noise covariance between each pair
% of units
% 

fnm = fullfile(saveDir, 'Fig5C.pdf');
doSave = false;

pts = [];
Cpairs = [];
c = 0;
for jj = 1:numel(dts)
    dtstr = dts{jj};
%     if ismember(dtstr, {'20150304a', '20150316c'})
%         continue;
%     end
    ix = ismember({pairs.dt}, dtstr);
    cpairs = pairs(ix);
    cpts = nan(numel(cpairs), 2);
    for ii = 1:numel(cpairs)
        cpair = cpairs(ii);
        cell1 = cells(ismember({cells.name}, cpair.cell1));
        cell2 = cells(ismember({cells.name}, cpair.cell2));
        cpair.cell1 = cell1; cpair.cell2 = cell2;
        rf1 = cell1.w; rf2 = cell2.w;
        cpts(ii,1) = rf1'*rf2;

        cpair.Ys_resAR = [cell1.YresAR cell2.YresAR];
        noiseCov = nancov(cpair.Ys_resAR);
        cpts(ii,2) = noiseCov(1,2);
        
%         cpts(ii,2) = cpts(ii,1);
%         cv = cov(rf1, rf2); cpts(ii,1) = cv(1,2);
%         cpts(ii,1) = cpair.rfCorr;
%         cpts(ii,1) = cpair.noiseCorrAR;
        
%         ixd = ~any(isnan(cpair.Ys_resAR),2);
%         assert(corr(cpair.Ys_resAR(ixd,1), cpair.Ys_resAR(ixd,2)) == ...
%             cpair.noiseCorrAR);
        Cpairs = [Cpairs cpair];
    end
    pts = [pts; cpts];        
end
Fig5C = plot.init;

xs = pts(:,1);
ys = pts(:,2);
% ixc = (xs >= prctile(xs, 1)) & (xs <= prctile(xs, 99));
% ixc = ixc & (ys >= prctile(ys, 1)) & (ys <= prctile(ys, 99));
% xs = xs(ixc); ys = ys(ixc);

set(gca, 'LineWidth', 2);
mdl = fitlm(xs, ys);
h = mdl.plot;
h(1).Marker = 'o';
h(1).Color = 'k';
h(1).MarkerSize = 5;
h(1).MarkerFaceColor = 0.5*ones(3,1);
h(1).LineWidth = 1.5;
h(2).LineWidth = 3;
h(2).Color = [0.8 0.2 0.2];
set(gca, 'TickDir', 'out');
legend off;
title('');
plot.setPrintSize(gcf, struct('width', 4, 'height', 3.4));

% ixp = abs(xs) < 3 & ys > 40;
% plot(xs(ixp), ys(ixp), 'ro');

% set(gca, 'LineWidth', 2);
% plot(pts(:,1), pts(:,2), 'ko', 'MarkerSize', 4.5, ...
%     'MarkerFaceColor', 0.5*ones(3,1));
xlabel(''); ylabel('');
xlabel('dot-product of subfields (w_1\cdotw_2)', ...
    'interpreter', 'tex', 'Color', 'k');
ylabel('noise covariance', 'Color', 'k');
% title(['r^2 = ' sprintf('%0.2f', corr(pts(:,1), pts(:,2)))]);
axis tight;

if doSave
    export_fig(Fig5C, fnm);
end

%% part 2. noise correlation predicts delta decoding (with groups of cells)

ks = 2:4; % 2 = pairs, 3 = triplets, 4 = quads, etc. (can include multiple)
fitToShuffled = true; % if false, fit to unshuffled data

ixCellsToKeep = figure.filterCellsAndPairs(cells, true, 0);
groups = tools.makeCellGroups(cells(ixCellsToKeep), ks);
dts = unique({groups.dt});

if fitToShuffled
    shufnm = 'shuffled';
else
    shufnm = 'unshuffled';
end
cnm = ['weights fit to ' shufnm];% '; k= ' num2str(ks)];

tic;
pts = [];
for jj = 1:numel(dts)
    dtstr = dts{jj};
    ix = ismember({groups.dt}, dtstr);
    cgroups = groups(ix);
    cpts = nan(numel(cgroups), 3);
    for ii = 1:numel(cgroups)
        cgroup = cgroups(ii);
        Ys_unshuf = cgroup.Ys_resAR;
        
        ixd = ~any(isnan(Ys_unshuf),2) & ~isnan(cgroup.stimdir);
        Ys_unshuf = Ys_unshuf(ixd,:);
        cgroup.stimdir = cgroup.stimdir(ixd,:);
        
        % define decoding weights based on aggregate tuning
        ws = -(2*(cgroup.targPrefs == 1)-1); ws = ws';
        
        grps = unique(cgroup.stimdir);
        Ys_shuf = nan(size(Ys_unshuf));
        for kk = 1:numel(grps)
            ixc = (cgroup.stimdir == grps(kk));
            Yc = Ys_unshuf(ixc,:);
            [~, idx] = sort(rand(size(Yc)),1);
            ytmp = nan(size(Yc));
            for ll = 1:size(Yc,2)
                ytmp(:,ll) = Yc(idx(:,ll),ll);
            end
            Ys_shuf(ixc,:) = ytmp;
        end
        
        Sigma = cov(Ys_unshuf);
        Sigma_diag = diag(diag(Sigma));
%         dPrime_shuf_hat = 1./(ws'*cov(Ys_shuf)*ws);
        dPrime_shuf_hat = 1./(ws'*Sigma_diag*ws);
        dPrime_raw_hat = 1./(ws'*Sigma*ws);
        cpts(ii,1) = sqrt(dPrime_raw_hat/dPrime_shuf_hat);
        
        if fitToShuffled
            Ysc = Ys_shuf;
        else
            Ysc = Ys_unshuf;
        end
        mdl = fitcdiscr(Ysc, cgroup.stimdir, 'DiscrimType', 'linear');
        K = mdl.Coeffs(1,2).Const; L = mdl.Coeffs(1,2).Linear;
        dPrime_shuf = tools.dprime(Ys_shuf*L + K, cgroup.stimdir);
        dPrime_raw = tools.dprime(Ys_unshuf*L + K, cgroup.stimdir);
        cpts(ii,2) = sqrt(dPrime_raw/dPrime_shuf);
        cpts(ii,3) = numel(ws); % whether pair, triplet, etc.
        
        cpts(ii,4) = dPrime_raw_hat;
        cpts(ii,5) = dPrime_shuf_hat;
        cpts(ii,6) = dPrime_raw;
        cpts(ii,7) = dPrime_shuf;
    end
    pts = [pts; cpts];
end
toc;

%%

pts = ptsFitToShuf;
fnm = fullfile(saveDir, 'Fig5D.pdf');
doSave = false;

Fig5D = plot.init;
set(gca, 'LineWidth', 2);
% plot(pts(:,1).^2, pts(:,2).^2, 'k.', 'MarkerSize', 20);
nks = unique(pts(:,3));
clrs = cbrewer('seq', 'Greens', 8);
clrs = clrs([4 6 8],:);
% clrs = clrs([3 5 8],:);

xlim([0.5 1.6]); ylim(xlim);
xlim([-1 1]); ylim(xlim);
plot(xlim, 0*[1 1], '--', 'Color', 0.8*ones(3,1), 'Linewidth', 2);
plot(0*[1 1], ylim, '--', 'Color', 0.8*ones(3,1), 'Linewidth', 2);

mdls = cell(numel(nks),1);
counts = nan(numel(nks),2,2);
for kk = numel(nks):-1:1
    ixc = pts(:,3) == nks(kk);
    clr = clrs(kk,:);
    xc = pts(ixc,1).^1;
    yc = pts(ixc,2).^1;
    
    xc = log(pts(ixc,4)) - log(pts(ixc,5));
    yc = log(pts(ixc,6)) - log(pts(ixc,7));
    
%     [bandwidth,density,xa,ya] = kde2d([xc yc]);
%     contour3(xa,ya,density,50);
%     surf(xa,ya,density,'LineStyle','none');
    
    scatter(xc, yc, 10, clr, 'o');%, 'filled');
%     mu = nanmean([xc yc]);
%     cpts = plot.gauss2dcirc([], 2, cov([xc yc]));
%     plot(mu(1)+cpts(1,:), mu(2)+cpts(2,:), '-', 'LineWidth', 2, ...
%         'Color', clr);
    
    counts(kk,1,1) = nanmean((pts(ixc,1) <= 1) & (pts(ixc,2) > 1));
    counts(kk,1,2) = nanmean((pts(ixc,1) > 1) & (pts(ixc,2) > 1));
    counts(kk,2,1) = nanmean((pts(ixc,1) <= 1) & (pts(ixc,2) <= 1));
    counts(kk,2,2) = nanmean((pts(ixc,1) > 1) & (pts(ixc,2) <= 1));
    round(100*squeeze(counts(kk,:,:)))
    100*(1 - sum(diag(squeeze(counts(kk,:,:))')))
    
    mdls{kk} = fitlm(pts(ixc,1), pts(ixc,2));
end
xlabel({'predicted log d'' ratio', ...
    '\leftarrow detrimental r_{sc}      beneficial r_{sc} \rightarrow'}, ...
    'Color', 'k');
ylabel({'actual log d'' ratio', ...
    '\leftarrow detrimental r_{sc}      beneficial r_{sc} \rightarrow'}, ...
    'Color', 'k');
ixc = ~any(isnan(pts(:,1:2)),2);
axis equal; axis tight;
set(gca, 'XTick', [-1 0 1]);
set(gca, 'YTick', get(gca, 'XTick'));
% title({cnm, ['r^2 = ' sprintf('%0.2f', corr(pts(ixc,1), pts(ixc,2)))]});
% xlim([0.5 1.6]); ylim(xlim);
if doSave
    export_fig(Fig5D, fnm);
end
