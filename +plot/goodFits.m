function vals = goodFits(fitstr, muCorrThresh, scoreThresh, scoreSdevThresh)
    if nargin < 4
        scoreSdevThresh = nan;
    end
    if nargin < 3
        scoreThresh = -0.2;
    end
    if nargin < 2
        muCorrThresh = 0.5;
    end
    dts = io.getDates('fits');
    c = 0;
    for ii = 1:numel(dts)
        dt = dts{ii};
        fs = io.loadFitsByDate(dt);
        nms = fieldnames(fs);
        for jj = 1:numel(nms)
            fit = fs.(nms{jj}).(fitstr){end};
            c = c+1;
            vals(c).type = nms{jj}(1:2);
            vals(c).score = mean(fit.scores);
%             vals(c).score = fit.score;
            mcf = fit.muCorrFolds;            
            vals(c).muCorr = min(mcf(abs(triu(mcf,1)-mcf)==0));
%             vals(c).scoreSdev = fit.scoreVarFolds; % scoreStdevFolds
            vals(c).scoreSdev = fit.scoreVarFolds*2/sqrt(numel(fit.scores));
            % mult by 2 above so that t/2 := score/2*scoreSdev > 1
            % since t scores greater than 2 are significant at p=0.05, right?
        end
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
    disp(['Ignoring ' num2str(sum(~inds)) ' fits.'])

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
