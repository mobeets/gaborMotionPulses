function vals = makeFitSummaries(dts, fitdir, fitstr)
    if nargin < 3
        fitstr = 'ASD';
    end
    if nargin < 2
        fitdir = 'fits';
    end
    if nargin < 1 || isnan(dts)
        dts = io.getDates(fitdir);
    end
    vals = struct([]);    
    for ii = 1:numel(dts)
        vals = [vals fitSummariesByDate(dts{ii}, fitdir, fitstr)];
    end
end

function vals = fitSummariesByDate(dt, fitdir, fitstr)
    d = io.loadDataByDate(dt); % can be replaced by fit.shape eventually
    fs = io.loadFitsByDate(dt, fitdir);
    nms = fieldnames(fs);
    
    vals = struct([]);
    for jj = 1:numel(nms)
        nm = strsplit(nms{jj}, '_');
        celltype = nm{1};
        fit = fs.(nms{jj}).(fitstr){end};

        val.dt = dt;
        val.type = celltype;
        val.name = [dt '-' nms{jj}];
        if numel(nm) > 1
            val.cellind = str2num(nm{2});
        else
            val.cellind = nan;
        end

        val.wf = fit.mu;            
        val.mu = reshape(fit.mu(1:end-1), d.ns, d.nt);        
        if ~strcmp(celltype, 'decision')
            val.ntrials = sum(~isnan(d.Y_all(:,jj)));
        else
            val.ntrials = sum(~isnan(d.R));
        end
        val.separability = getSeparability(val.mu);
        val.score_mean = mean(fit.scores);
        val.score = fit.score;
        mcf = fit.muCorrFolds;
        val.muCorr = min(mcf(abs(triu(mcf,1)-mcf)==0));
        if isfield(fit, 'scoreVarFolds')
            val.scoreSdev = fit.scoreVarFolds; % scoreStdevFolds
        else
            val.scoreSdev = fit.scoreStdvFolds; % scoreStdevFolds
        end
        val.tScoreDenom = val.scoreSdev*2/sqrt(numel(fit.scores));
        % p < 0.05 if val.score_mean / val.tScoreDenom > 1
        
        vals = [vals val];
    end
    vals = correlateCellAndDecisionMu(vals);
end

function val = getSeparability(mu)    
    [~,s,~] = svd(mu);
    s = diag(s);
    val = s(1)/sum(s);
end

function vals = correlateCellAndDecisionMu(vals)
    decinds = strcmp({vals.type}, 'decision');
    dec = vals(decinds);
    for ii = 1:numel(vals)
        r = corrcoef(dec.mu, vals(ii).mu);
        vals(ii).decisionCorrelation = r(2);
    end
end
