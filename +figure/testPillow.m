%% load cell pairs

pairs_all = tools.makeCellPairs(cells);
ixPairsToKeep = figure.filterCellsAndPairs(pairs_all, false);
pairs = pairs_all(ixPairsToKeep);
dts = unique({pairs.dt});

%% part 1. RFs predict noise correlation

% current status: I think this is done. we don't need to do triplets or
% more (i.e., pairs is fine) because the claim is only that the dot-product
% of the RF maps is proportional to the noise covariance between each pair
% of units
% 

pts = [];
for jj = 1:numel(dts)
    dtstr = dts{jj};
    ix = ismember({pairs.dt}, dtstr);
    cpairs = pairs(ix);
    cpts = nan(numel(cpairs), 2);
    for ii = 1:numel(cpairs)
        cpair = cpairs(ii);
        cell1 = cells(ismember({cells.name}, cpair.cell1));
        cell2 = cells(ismember({cells.name}, cpair.cell2));
        rf1 = cell1.w; rf2 = cell2.w;
        cpts(ii,1) = rf1'*rf2;        

        noiseCov = nancov(cpair.Ys_resAR);
        cpts(ii,2) = noiseCov(1,2);
        
%         cpts(ii,2) = cpts(ii,1);
%         cv = cov(rf1, rf2); cpts(ii,1) = cv(1,2);
%         cpts(ii,1) = cpair.rfCorr;
%         cpts(ii,1) = cpair.noiseCorrAR;
        
        ixd = ~any(isnan(cpair.Ys_resAR),2);
        assert(corr(cpair.Ys_resAR(ixd,1), cpair.Ys_resAR(ixd,2)) == ...
            cpair.noiseCorrAR);
    end
    pts = [pts; cpts];    
end
plot.init;
plot(pts(:,1), pts(:,2), 'k.', 'MarkerSize', 20);
xlabel('predicted noise cov. (using STRFs)');
ylabel('actual noise cov. (using residuals)');
title(['r^2 = ' sprintf('%0.2f', corr(pts(:,1), pts(:,2)))]);

%% part 2. noise correlation predicts delta decoding (with groups of cells)

ks = 2:4; % 2 = pairs, 3 = triplets, 4 = quads, etc. (can include multiple)
fitToShuffled = true; % if false, fit to unshuffled data

ixCellsToKeep = figure.filterCellsAndPairs(cells, false);
cellGroups = tools.makeCellGroups(cells(ixCellsToKeep), ks);
ixGroupsToKeep = figure.filterCellsAndPairs(cellGroups, false);
groups = cellGroups(ixGroupsToKeep);
dts = unique({groups.dt});

if fitToShuffled
    shufnm = 'shuffled';
else
    shufnm = 'unshuffled';
end
cnm = ['weights fit to ' shufnm];% '; k= ' num2str(ks)];
% todo next:
% - ideally I think we'd like the prediction to be in terms of
% delta-percent correct, since that's what we use elsewhere in the paper
% - need to do cross-val on model fitting

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
        
        % n.b. need to do crossval
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
    end
    pts = [pts; cpts];
end

%%

plot.init;
% plot(pts(:,1).^2, pts(:,2).^2, 'k.', 'MarkerSize', 20);
nks = unique(pts(:,3));
clrs = [0.8 0.2 0.2; 0.2 0.2 0.8; 0.2 0.8 0.2];
for kk = numel(nks):-1:1
    ixc = pts(:,3) == nks(kk);
    clr = clrs(kk,:);
    plot(pts(ixc,1).^1, pts(ixc,2).^1, '.', 'Color', clr, ...
        'MarkerSize', 1);
end
xlabel({'predicted d-prime ratio', ...
    '\leftarrow r_{sc} hurts         r_{sc} helps \rightarrow'});
ylabel({'actual d-prime ratio', ...
    '\leftarrow r_{sc} hurts         r_{sc} helps \rightarrow'});
ixc = ~any(isnan(pts(:,1:2)),2);
axis equal;
xlim([0.5 1.6]); ylim(xlim);
set(gca, 'XTick', [0.5 1 1.5]); set(gca, 'YTick', [0.5 1 1.5]);
title({cnm, ['r^2 = ' sprintf('%0.2f', corr(pts(ixc,1), pts(ixc,2)))]});
