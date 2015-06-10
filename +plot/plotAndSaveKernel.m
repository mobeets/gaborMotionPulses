function plotAndSaveKernel(obj, data, figdir)
% 
% plots the kernel and saves to png
% 
    fig_fnfcn = @(tag, ext) fullfile(figdir, [tag '.' ext]);
    fig_svfcn = @(fig, tag, ext) hgexport(fig, fig_fnfcn(tag, ext), ...
        hgexport('factorystyle'), 'Format', ext);
    fig_lblfcn = @(lbl, sc) [lbl ' sc=' num2str(sprintf('%.2f', sc))];
    
    if ~isempty(figdir)
        wf = obj.mu;
        if prod(obj.shape) < wf
            wf = wf(1:end-1);
        end
        fig = plot.plotKernel(data.Xxy, reshape(wf, data.ns, data.nt), ...
            nan, fig_lblfcn(obj.label, obj.score));
        fig_svfcn(fig, obj.label, 'png');        
    end
end
