function fits = fitSTRF(data, foldinds, mapFcn, mlFcn, bmapFcn, scFcn, ...
    lbs, ubs, ns, figdir, lbl, ifold, mask)
% fits = fitSTRF(data, foldinds, mapFcn, mlFcn, scFcn, ...
%     lbs, ubs, ns, figdir, lbl, ifold, mask)
% 
% wrapper for fitting and plotting space-time receptive fields
%   using either gridding of hyperparameter space, or a zooming grid search
% 
% - data.X, data.Y - regressors and response
% - mapFcn, mlFcn, bmapFcn - for ASD, ML, and bilinear ASD, respectively
% - foldinds, scFcn, lbs, ubs, ns - for cross-validation and hyperparams
% - figdir, lbl, ifold - for plotting
% - mask - which fitting/search methods to use
% 

    fitMap = mask(1);
    fitMapGridSearch = mask(2);
    fitML = mask(3);
    fitMapBilinear = mask(4);

    isLog = true; % lbs and ubs are in logspace    

    % MAP estimate on each hyper in hypergrid
    if fitMap
        hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
        obj = reg.cvMaxScoreGrid(data.X, data.Y, hypergrid, mapFcn, {}, ...
            scFcn, {}, foldinds);
        fits.ASD = addFigure(obj, data, [lbl '-ASD'], figdir, ifold);
    end

    % MAP estimate on hypergrid with grid search
    if fitMapGridSearch
        obj = reg.cvMaxScoreGridSearch(data.X, data.Y, lbs, ubs, ns, ...
            mapFcn, {}, scFcn, {}, foldinds, isLog);
        fits.ASD_gs = addFigure(obj, data, [lbl '-ASD-gs'], figdir, ifold);
    end

    % ML estimate on each hyper in hypergrid
    if fitML
        obj = reg.cvMaxScoreGrid(data.X, data.Y, [nan nan nan], mlFcn, ...
            {}, scFcn, {}, foldinds);
        fits.ML = addFigure(obj, data, [lbl '-ML'], figdir, ifold);
    end
    
    % MAP estimate using separable space and time weights
    if fitMapBilinear
        hypergrid = exp(tools.gridCartesianProduct(lbs(1:end-1), ...
            ubs(1:end-1), ns(1:end-1))); % ignore hyper for time smoothing
        obj = reg.cvMaxScoreGrid(data.X, data.Y, hypergrid, bmapFcn, ...
            {}, scFcn, {}, foldinds);
        fits.ASD_b = addFigure(obj, data, [lbl '-ASDb'], figdir, ifold);
    end

end
%%

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
    if ~isempty(figdir)
        wf = obj.mus{ifold}; sc = obj.scores(ifold);
        fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), ...
            nan, fig_lblfcn(obj.label, 1, sc));
        fig_svfcn(fig, obj.label, 'png');        
    end
end
