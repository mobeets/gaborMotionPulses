function plotAndSaveSummaries(fitdir, outdir, isNancy, figext)
% plotAndSaveSummaries(fitdir, outdir, figext)
% 
% generates summary plots for all subfolders in fitdir
% saves figures to outdir with filetype specified by figext
% 
    if nargin < 4
        figext = 'png';
    end
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    else
        warning(['About to overwrite folder "' outdir '"']);
        x = input('Continue? ', 's');
        if ~strcmpi(x(1), 'y')
            return;
        end
    end

    dts = io.getDates(fitdir);
    for ii = 1:numel(dts)
        dt = dts{ii};
        disp(dt);
        curdir = fullfile(outdir, dt);
        if ~exist(curdir, 'dir')
            mkdir(curdir);
        end        
        fig = plot.quickPmfByDate(dt, isNancy);
        plot.saveFig(fig, 'pmf', curdir, figext);
        fig = plot.kernelsByDate(dt, fitdir, isNancy);
        plot.saveFig(fig, 'kernels', curdir, figext);
        plot.summaryByCell(dt, nan, isNancy, fitdir, curdir, figext);
        close all
    end
    
%     fig = plot.scoreVsHyperparam(1, fitdir);
    fig = plot.separabilityVsScore(fitdir);
    plot.saveFig(fig, 'sepVsScore', outdir, figext);
    close all
end
