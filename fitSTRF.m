function fits = fitSTRF(data, ML, MAP, BMAP, eMAP, scoreFcn, hyperOpts, ...
    figdir, lbl, mask, foldinds, evalinds)
% fits = fitSTRF(data, mapFcn, mlFcn, scFcn, ...
%     lbs, ubs, ns, figdir, lbl, mask, foldinds, evalinds)
% 
% wrapper for fitting and plotting space-time receptive fields
%   using either gridding of hyperparameter space, or a zooming grid search
% 
% - data.X, data.Y - regressors and response
% - mapFcn, mlFcn, bmapFcn - for ASD, ML, and bilinear ASD, respectively
% - foldinds, scFcn, lbs, ubs, ns - for cross-validation and hyperparams
% - figdir, lbl - for plotting
% - mask - which fitting/search methods to use
% - evalinds - trials to use for finding hyperparams
% 

    fitMap = mask(1);
    fitMapGridSearch = mask(2);
    fitML = mask(3);
    fitMapBilinear = mask(4);
    fitMapEviOpt = mask(5);
        
    [X, Y, foldinds, evalinds] = io.dropTrialsIfYIsNan(data.X, data.Y, ...
        foldinds, evalinds);

    % MAP estimate on each hyper in hypergrid
    if fitMap
        % hypergrid, exp on appropriate parts as dictated by isLog
        lbs = hyperOpts.lbs; ubs = hyperOpts.ubs; ns = hyperOpts.ns;
        hypergrid = tools.gridCartesianProduct(lbs, ubs, ns);
        isLog = repmat(hyperOpts.isLog, size(hypergrid,1), 1);
        hypergrid = isLog.*(exp(hypergrid)) + (1-isLog).*hypergrid;
        
        obj = reg.cvMaxScoreGrid(X, Y, MAP, scoreFcn, hypergrid, ...
            foldinds, evalinds, 'grid');
        fits.ASD = plot.plotAndSaveKernel(obj, data, [lbl '-ASD-g'], figdir);
    end

    % MAP estimate on hypergrid with grid search
    if fitMapGridSearch
        obj = reg.cvMaxScoreGrid(X, Y, MAP, scoreFcn, nan, ...
            foldinds, evalinds, 'grid-search', hyperOpts);
        fits.ASD_gs = plot.plotAndSaveKernel(obj, data, [lbl '-ASD-gs'], figdir);
    end

    % ML estimate
    if fitML
        obj = reg.fitToEvalinds(X, Y, evalinds, ML, scoreFcn, nan);
        obj.score_noCV = obj.score;
        objCV = reg.cvMaxScoreGrid(X, Y, ML, scoreFcn, [nan nan], ...
            foldinds, evalinds, 'grid');
        obj.score = objCV.score; obj.scores = objCV.scores;
        obj.score_cvMean = mean(obj.scores);
        obj.score_cvStd = std(obj.scores);
        fits.ML = plot.plotAndSaveKernel(obj, data, [lbl '-ML'], figdir);
    end
    
    % MAP estimate using separable space and time weights
    if fitMapBilinear
        % fit hyperparameter under fully smooth temporal filter
        sMAP = @(hyper) asd.fitHandle(hyper, data.D, ...
            'gauss', 'evi', struct('fullTemporalSmoothing', true));
        obj = reg.fitToEvalinds(X, Y, evalinds, sMAP, scoreFcn, nan);
        hyper = obj.hyper(1:end-1);

        % fit
        obj = reg.fitToEvalinds(X, Y, evalinds, BMAP, scoreFcn, hyper);
        obj.score_noCV = obj.score;
        
        % cv to score
        objCV = reg.cvMaxScoreGrid(X, Y, BMAP, scoreFcn, hyper', ...
            foldinds, evalinds, 'grid');
        obj.score = objCV.score; obj.scores = objCV.scores;
        obj.score_cvMean = mean(obj.scores);
        obj.score_cvStd = std(obj.scores);
        fits.ASD_b = plot.plotAndSaveKernel(obj, data, [lbl '-ASDb'], figdir);
    end
    
    % MAP estimate using evidence optimization (requires gaussian ll)
    if fitMapEviOpt
        obj = reg.fitToEvalinds(X, Y, evalinds, eMAP, scoreFcn, nan);
        obj.score_noCV = obj.score;
        objCV = reg.cvMaxScoreGrid(X, Y, MAP, scoreFcn, obj.hyper', ...
            foldinds, evalinds, 'grid');
        obj.score = objCV.score; obj.scores = objCV.scores;
        obj.score_cvMean = mean(obj.scores);
        obj.score_cvStd = std(obj.scores);
        fits.ASD = plot.plotAndSaveKernel(obj, data, [lbl '-ASD'], figdir);
    end

end
