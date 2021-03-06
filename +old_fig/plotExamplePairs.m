% doSave = false;
% fnm = fullfile('data', 'figs', 'summary', [fitnm '_examples.pdf']);

scs = pairs;
% minDprime = 0.4;
% ix = ([scs.cell1_dPrime] > minDprime) & ([scs.cell2_dPrime] > minDprime);
ix = [scs.minAbsDprime] > 0.4;
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
ixC = samePool & posNoiseCorr & [scs.scoreGainWithCorrs] < -0.03;
psC = scs(ixC);
psC = psC(4);

% same pools, helpful noise corrs
ixD = samePool & posNoiseCorr & [scs.scoreGainWithCorrs] > 0.02;
psD = scs(ixD);

%%

doSave = true;
saveDir = 'data/figs/summary/cell-pairs';
if doSave && ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

ps = [psC psB psA];
nms = {'same_pool_pos_noisecorr', 'diff_pool_neg_noisecorr', ...
    'diff_pool_pos_noisecorr_1', 'diff_pool_pos_noisecorr_2'};

% ps = psA;
% nms = {'diff_pool_pos_noisecorr'};
for ii = 1:numel(ps)
    plot.init;
    
    % show scatter
    if ~doSave
        subplot(2,2,1); hold on;
    end
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
    if doSave
        plot.setPrintSize(gcf, struct('width', 3.5, 'height', 3));
        cnm = [nms{ii} '_scatter'];
        fnm = fullfile(saveDir, [cnm '.pdf']);
        export_fig(gcf, fnm);
    end
    
    cell1 = cells(strcmpi({cells.name}, ps(ii).cell1));
    cell2 = cells(strcmpi({cells.name}, ps(ii).cell2));
    sz = 120;
    
    % show RF of y-axis cell
    if ~doSave
        subplot(2,2,2); hold on;
    else
        plot.init;
    end
    plot.rfSingle(cell2.Xxy, cell2.wsep.spatial_RF, nan, sz);
    if doSave
        plot.setPrintSize(gcf, struct('width', 3, 'height', 3));
        cnm = [nms{ii} '_y_' cell2.name];
        fnm = fullfile(saveDir, [cnm '.pdf']);
        export_fig(gcf, fnm);
    else
        title(cell2.name);
    end
    
    % show RF of x-axis cell
    if ~doSave
        subplot(2,2,3); hold on;
    else
        plot.init;
    end
    plot.rfSingle(cell1.Xxy, cell1.wsep.spatial_RF, nan, sz);
    if doSave
        plot.setPrintSize(gcf, struct('width', 3, 'height', 3));
        cnm = [nms{ii} '_x_' cell1.name];
        fnm = fullfile(saveDir, [cnm '.pdf']);
        export_fig(gcf, fnm);
    else
        title(cell1.name);
    end
    
    if ~doSave
        plot.setPrintSize(gcf, struct('width', 7, 'height', 6));
    end

end
