function saveFigure(name, basedir, fig, ext)
    if nargin < 4
        ext = 'png';
    end
    if nargin < 3
        fig = gcf;
    end
    if nargin < 2
        basedir = 'tmp';
    end
    if ~exist(basedir, 'dir')
        mkdir(basedir);
    end
    fn = fullfile(basedir, [name '.' ext]);
    hgexport(fig, fn, hgexport('factorystyle'), 'Format', ext);
end
