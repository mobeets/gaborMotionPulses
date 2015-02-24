function fits = innermain(data, foldinds, mapFcn, mlFcn, scFcn, lbs, ubs, ns, figdir, lbl, ifold, mask)
    doASD = mask(1);
    doASD_gs = mask(2);
    doML = mask(3);

    isLog = true; % lbs and ubs are in logspace
    hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
    fig_fnfcn = @(tag, ext) fullfile(figdir, [tag '.' ext]);
    fig_svfcn = @(fig, tag, ext) hgexport(fig, fig_fnfcn(tag, ext), hgexport('factorystyle'), 'Format', ext);
    fig_lblfcn = @(lbl, ifold, sc) [lbl ' f' num2str(ifold) ' sc=' num2str(sprintf('%.2f', sc))];
    wf_fcn = @(wf, ns, nt) reshape(wf(1:end-1), ns, nt);

    % ASD on hypergrid
    if doASD
        obj = reg.cvMaxScoreGrid(data, hypergrid, mapFcn, {}, ...
            scFcn, {}, foldinds);
        obj.label = [lbl '-ASD'];
        if ~isempty(figdir)
            wf = obj.mus{ifold}; sc = obj.scores(ifold);
            fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), nan, ...
                fig_lblfcn(obj.label, 1, sc));
            fig_svfcn(fig, obj.label, 'png');
        end
        fits.ASD = obj;
    end

    % ASD on hypergrid with grid search
    if doASD_gs
        obj = reg.cvMaxScoreGridSearch(data, lbs, ubs, ns, mapFcn, {}, ...
            scFcn, {}, foldinds, isLog);
        obj.label = [lbl '-ASD-gs'];
        if ~isempty(figdir)
            wf = obj.mus{ifold}; sc = obj.scores{ifold};
            fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), nan, ...
                fig_lblfcn(obj.label, 1, sc));
            fig_svfcn(fig, obj.label, 'png');
        end
        fits.ASD_gs = obj;
    end

    % ML
    if doML
        obj = reg.cvMaxScoreGrid(data, [nan nan nan], mlFcn, {}, ...
            scFcn, {}, foldinds);
        obj.label = [lbl '-ML'];
        if ~isempty(figdir)
            wf = obj.mus{ifold}; sc = obj.scores(ifold);
            fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), nan, ...
                fig_lblfcn(obj.label, 1, sc));
            fig_svfcn(fig, obj.label, 'png');        
        end
        fits.ML = obj;
    end

end
