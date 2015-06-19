function plotAndSaveOverlays(fitdir, isNancy, figdir, cellType, fitstr, ...
    cellsPerFit)

    if nargin < 5
        fitstr = 'ASD';
    end
    if nargin < 6
        cellsPerFit = 1;
    end
    
    dts = io.getDates(fitdir);
    for ii = 1:numel(dts)
        dt = dts{ii};
        d = io.loadDataByDate(dt, isNancy);        
        if isempty(fieldnames(d))
            continue;
        end
        inds = arrayfun(@(n) n.dPrime > 0.5 & ...
            strcmp(n.brainArea, cellType), [d.neurons{:}]);
        if strcmp(cellType, 'LIP')
            inds2 = arrayfun(@(n) ~isempty(n.delayedsaccades), [d.neurons{:}]);
        elseif strcmp(cellType, 'MT')
            inds2 = arrayfun(@(n) ~isempty(n.hyperflow), [d.neurons{:}]);
        elseif strcmp(cellType, 'decision')
            fs = io.loadFitsByDate(dt, fitdir);
            f = fs.decision.(fitstr);
            figure;
            plot.plotSaccadeKernelOverlay(d.stim, nan, f{end}, true, true);
            plot.saveFigure([dt '-decision'], figdir);
            continue;
        end
        ix = 1:numel(d.neurons);
        cellinds = ix(inds & inds2);
        nblks = floor(numel(cellinds)/cellsPerFit)+1;
        for jj = 1:nblks
            from = (jj-1)*cellsPerFit+1;
            to = ((jj-1)*cellsPerFit)+cellsPerFit;
            cix = cellinds(from:min(to,end));
            if isempty(cix)
                continue;
            end
            fig = plot.plotAllSaccadeKernelOverlays(dt, fitdir, isNancy, ...
                fitstr, cix, 2);
            set(fig, 'Position', [100, 100, 1500, 350]);
            lbl = [dt '-' cellType '-' strjoin({num2str(cix)}, ',')];
            plot.saveFigure(lbl, figdir);
        end
    end
    close all
end
