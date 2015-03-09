function saveFig(fig, name, outdir, ext)
% function saveFig(fig, nm, outdir, ext)
% 
% fig - figure handle
% name - filename for saved figure, without extension
% outdir - output directory for saved figure
% ext - extension of figure, e.g., 'png'
% 
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, fullfile(outdir, name), ['-d' ext], '-r0');
end
