function vals = makeFitSummaries(fitdir, isNancy, fitstr, dts)      
    if nargin < 1
        fitdir = 'fits';
    end
    if nargin < 2
        isNancy = false;
    end
    if nargin < 3
        fitstr = 'ASD';
    end
    if nargin < 4 || isnan(dts)
        dts = io.getDates(fitdir);
    end    
    
    vals = struct([]);    
    for ii = 1:numel(dts)
        vals = [vals fitSummariesByDate(fitdir, dts{ii}, fitstr, isNancy)];
    end
end

function vals = fitSummariesByDate(fitdir, dt, fitstr, isNancy)
    d = io.loadDataByDate(dt, isNancy);
    fs = io.loadFitsByDate(dt, fitdir);
    nms = fieldnames(fs);
    
    vals = struct([]);
    for jj = 1:numel(nms)
        nm = strsplit(nms{jj}, '_');
        celltype = nm{1};
        if ~isfield(fs.(nms{jj}), fitstr)
            continue;
        end
        fit = fs.(nms{jj}).(fitstr){end};

        val.dt = dt;
        val.type = celltype;
        val.name = [dt '-' nms{jj}];        
        if numel(nm) > 1
            val.cellind = str2num(nm{2});
            cip = ['-' nm{2}];
        else
            val.cellind = nan;
            cip = '';
        end
        
        val.isNancy = isNancy;
        val.pngname = fullfile(val.dt, [val.type cip '-ASD.png']);
        
        val.wf = fit.mu;
        val.mu = reshape(fit.mu(1:end-1), d.ns, d.nt);

        if isfield(fit, 'hyper') && all(isnan(fit.hyper))
            fit.hyper = nan(3,1);
        end
        val.hyper_ro = fit.hyper(1);
        val.hyper_delta_space = fit.hyper(end-1);
        val.hyper_delta_time = fit.hyper(end);
        
        if ~strcmp(celltype, 'decision')
            val.ntrials = sum(~isnan(d.Y_all(:,val.cellind)));
            val.dPrime = d.neurons{val.cellind}.dPrime;
            val.hyper_ssq = fit.hyper(2);
        else
            val.dPrime = nan;
            val.ntrials = sum(~isnan(d.R));
            val.hyper_ssq = nan;
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
