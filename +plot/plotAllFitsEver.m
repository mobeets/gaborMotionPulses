function plotAllFitsEver(basedir, isNancy, dts, fitType)
    if nargin < 4
        fitType = 'ASD';
    end
    fitdir = fullfile(basedir, 'fits');
    figdir = fullfile(basedir, 'figs');
    outdir = fullfile(figdir, 'all');
    if nargin < 3
        dts = io.getDates(fitdir);
    end
    if nargin < 2
        isNancy = true;
    end
    
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    for kk = 1:numel(dts)
        dt = dts{kk};
        fs = io.loadFitsByDate(dt, fitdir);
        d = io.loadDataByDate(dt, isNancy);

        fns = fieldnames(fs);
        for ii = 1:numel(fns)
            fn = fns{ii};
            fs0 = fs.(fn).(fitType);
            for jj = 1:numel(fs0)
                f = fs0{jj};
                lbl = [dt ' - ' f.label ' - #' num2str(jj) ' - ' ...
                    sprintf('%0.0f, ', f.hyper)];
                plot.plotAndSaveKernel(f, d, lbl, outdir);
            end
            close all;
        end        
    end
end
