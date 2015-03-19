function vals = fitMetrics(vals, muCorrThresh, scoreThresh, scoreSdevThresh)
    if nargin < 4
        scoreSdevThresh = nan;
    end
    if nargin < 3
        scoreThresh = nan;
    end
    if nargin < 2
        muCorrThresh = 0.5;
    end
    
    inds = true(numel(vals),1)';
    if ~isnan(scoreSdevThresh)
        inds = inds & ([vals.scoreVar] >= scoreSdevThresh);
    end
    if ~isnan(scoreThresh)
        inds = inds & ([vals.score] >= scoreThresh);
    end
    if ~isnan(muCorrThresh)
        inds = inds & ([vals.muCorr] >= muCorrThresh);
    end
    vals = vals(inds);
    disp(['Ignoring ' num2str(sum(~inds)) ' fits.']);

    figure; hold on;
    xs = [vals.score];
    ys = [vals.muCorr];
    zs = [vals.scoreSdev];
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    for ii = 1:numel(categs)
        ind = strcmp(categs(ii), {vals.type});
        scatter3(xs(ind), ys(ind), zs(ind), 40, clrs(ii,:), 'filled', ...
            'DisplayName', categs(ii));
    end
    xlabel('score'); ylabel('mu-corr'); zlabel('score-sdev');
    legend(categs, 'Location', 'NorthEastOutside');
end
