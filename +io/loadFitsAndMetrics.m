function vals = loadFitsAndMetrics(fitstr, ignoredCellTypes)
    if nargin < 2
        ignoredCellTypes = {};
    end
    c = 0;
    dts = io.getDates('fits');    
    for ii = 1:numel(dts)
        dt = dts{ii};
        d = io.loadDataByDate(dt); % can be replaced by fit.shape eventually
        fs = io.loadFitsByDate(dt);
        nms = fieldnames(fs);
        for jj = 1:numel(nms)
            nm = strsplit(nms{jj}, '_');
            celltype = nm{1};
            if ismember(celltype, ignoredCellTypes)
                continue;
            end
            fit = fs.(nms{jj}).(fitstr){end};
            
            c = c+1;
            vals(c).dt = dt;
            vals(c).type = celltype;
            vals(c).name = [dt '-' nms{jj}];
            
            vals(c).wf = fit.mu;
            vals(c).mu = reshape(fit.mu(1:end-1), d.ns, d.nt);
            vals(c).score_mean = mean(fit.scores);
            vals(c).score = fit.score;
            mcf = fit.muCorrFolds;
            vals(c).muCorr = min(mcf(abs(triu(mcf,1)-mcf)==0));
            vals(c).scoreSdev = fit.scoreVarFolds; % scoreStdevFolds
            vals(c).tScoreDenom = fit.scoreVarFolds*2/sqrt(numel(fit.scores));
            % p < 0.05 if vals(c).score_mean / vals(c).tScoreDenom > 1
        end
    end
end
