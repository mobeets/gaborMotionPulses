function vals = fitSummariesByDate(dt, fitstr)
    d = io.loadDataByDate(dt); % can be replaced by fit.shape eventually
    fs = io.loadFitsByDate(dt);
    nms = fieldnames(fs);
    
    vals = struct([]);
    for jj = 1:numel(nms)
        nm = strsplit(nms{jj}, '_');
        celltype = nm{1};
        fit = fs.(nms{jj}).(fitstr){end};

        val.dt = dt;
        val.type = celltype;
        val.name = [dt '-' nms{jj}];

        val.wf = fit.mu;            
        val.mu = reshape(fit.mu(1:end-1), d.ns, d.nt);
        val.separability = getSeparability(val.mu);
        val.score_mean = mean(fit.scores);
        val.score = fit.score;
        mcf = fit.muCorrFolds;
        val.muCorr = min(mcf(abs(triu(mcf,1)-mcf)==0));
        val.scoreSdev = fit.scoreVarFolds; % scoreStdevFolds
        val.tScoreDenom = fit.scoreVarFolds*2/sqrt(numel(fit.scores));
        % p < 0.05 if val.score_mean / val.tScoreDenom > 1
        
        vals = [vals val];
    end
end

function val = getSeparability(mu)    
    [~,s,~] = svd(mu);
    s = diag(s);
    val = s(1)/sum(s);
end
