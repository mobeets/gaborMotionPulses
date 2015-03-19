function fits = fitSTRF(data, ML, MAP, BMAP, scoreFcn, hyperOpts, ...
    figdir, lbl, ifold, mask, foldinds, evalinds)
% fits = fitSTRF(data, mapFcn, mlFcn, scFcn, ...
%     lbs, ubs, ns, figdir, lbl, ifold, mask, foldinds, evalinds)
% 
% wrapper for fitting and plotting space-time receptive fields
%   using either gridding of hyperparameter space, or a zooming grid search
% 
% - data.X, data.Y - regressors and response
% - mapFcn, mlFcn, bmapFcn - for ASD, ML, and bilinear ASD, respectively
% - foldinds, scFcn, lbs, ubs, ns - for cross-validation and hyperparams
% - figdir, lbl, ifold - for plotting
% - mask - which fitting/search methods to use
% - evalinds - trials to use for finding hyperparams
% 

    fitMap = mask(1);
    fitMapGridSearch = mask(2);
    fitML = mask(3);
    fitMapBilinear = mask(4);
    
    lbs = hyperOpts.lbs; ubs = hyperOpts.ubs; ns = hyperOpts.ns;
    [X, Y, foldinds, evalinds] = dropTrialsIfYIsNan(data.X, data.Y, ...
        foldinds, evalinds);

    % MAP estimate on each hyper in hypergrid
    if fitMap
        hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
        obj = reg.cvMaxScoreGrid(X, Y, MAP, scoreFcn, hypergrid, ...
            foldinds, evalinds, 'grid');
        fits.ASD = addFigure(obj, data, [lbl '-ASD'], figdir, ifold);
    end

    % MAP estimate on hypergrid with grid search
    if fitMapGridSearch
        obj = reg.cvMaxScoreGrid(X, Y, MAP, scoreFcn, nan, ...
            foldinds, evalinds, 'grid-search', hyperOpts);
        fits.ASD_gs = addFigure(obj, data, [lbl '-ASD-gs'], figdir, ifold);
    end

    % ML estimate on each hyper in hypergrid
    if fitML
        obj = reg.cvMaxScoreGrid(X, Y, ML, scoreFcn, [nan nan nan], ...
            foldinds, evalinds, 'grid');
        fits.ML = addFigure(obj, data, [lbl '-ML'], figdir, ifold);
    end
    
    % MAP estimate using separable space and time weights
    if fitMapBilinear
        hypergrid = exp(tools.gridCartesianProduct(lbs(1:end-1), ...
            ubs(1:end-1), ns(1:end-1))); % ignore hyper for time smoothing
        obj = reg.cvMaxScoreGrid(X, Y, BMAP, scoreFcn, hypergrid, ...
            foldinds, evalinds, 'grid');
        fits.ASD_b = addFigure(obj, data, [lbl '-ASDb'], figdir, ifold);
    end

end
%%

function [X, Y, foldinds, evalinds] = dropTrialsIfYIsNan(X, Y, ...
    foldinds, evalinds)
    inds = isnan(Y);
    X = X(~inds,:,:);
    Y = Y(~inds,:);
    
    % remove foldinds relative to longer inds list
    if numel(inds) > numel(foldinds)
        f_inds = nan(numel(inds),1);
        f_inds(evalinds) = foldinds;
        f_inds(inds) = NaN;
        foldinds = f_inds(~isnan(f_inds));
    else
        foldinds = foldinds(~inds);
    end
    
    evalinds = evalinds(~inds);
end

function obj = addFigure(obj, data, lbl, figdir, ifold)
% 
% plots the kernel and saves to png
% 
    fig_fnfcn = @(tag, ext) fullfile(figdir, [tag '.' ext]);
    fig_svfcn = @(fig, tag, ext) hgexport(fig, fig_fnfcn(tag, ext), ...
        hgexport('factorystyle'), 'Format', ext);
    fig_lblfcn = @(lbl, ifold, sc) [lbl ' f' num2str(ifold) ...
        ' sc=' num2str(sprintf('%.2f', sc))];
    wf_fcn = @(wf, ns, nt) reshape(wf(1:end-1), ns, nt);
    
    obj.label = lbl;
    obj.shape = [data.ns data.nt];
    if ~isempty(figdir)
        wf = obj.mus{ifold}; sc = obj.scores(ifold);
        fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), ...
            nan, fig_lblfcn(obj.label, 1, sc));
        fig_svfcn(fig, obj.label, 'png');        
    end
end
