function obj = plotAndSaveKernel(obj, data, lbl, figdir)
% 
% plots the kernel and saves to png
% 
    fig_fnfcn = @(tag, ext) fullfile(figdir, [tag '.' ext]);
    fig_svfcn = @(fig, tag, ext) hgexport(fig, fig_fnfcn(tag, ext), ...
        hgexport('factorystyle'), 'Format', ext);
    fig_lblfcn = @(lbl, sc) [lbl ' sc=' num2str(sprintf('%.2f', sc))];
    wf_fcn = @(wf, ns, nt) reshape(wf(1:end-1), ns, nt);
    
    obj.label = lbl;
    obj.shape = [data.ns data.nt];
    if ~isempty(figdir)
        wf = obj.mu;
        sc = obj.score;
        fig = plot.plotKernel(data.Xxy, wf_fcn(wf, data.ns, data.nt), ...
            nan, fig_lblfcn(obj.label, sc));
        fig_svfcn(fig, obj.label, 'png');        
    end
end
