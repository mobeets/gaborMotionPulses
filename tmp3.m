
isNancy = true;
fitdir = 'data/cur-nancy/fits';
fitstr = 'ASD';
cellType = 'LIP';

figdir = 'figures/LIP';
fig_fnfcn = @(tag) fullfile(figdir, [tag '-plain.png']);
fig_svfcn = @(fig, tag) hgexport(fig, fig_fnfcn(tag), ...
    hgexport('factorystyle'), 'Format', 'png');

cellsPerFit = 1;

dts = io.getDates(fitdir);
for ii = 1:numel(dts)
    dt = dts{ii};
    d = io.loadDataByDate(dt, isNancy);
    inds = arrayfun(@(n) n.dPrime > 0.5 & ...
        strcmp(n.brainArea, cellType), [d.neurons{:}]);
    if strcmp(cellType, 'LIP')
        inds2 = arrayfun(@(n) ~isempty(n.delayedsaccades), [d.neurons{:}]);
    else
        inds2 = arrayfun(@(n) ~isempty(n.hyperflow), [d.neurons{:}]);
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
            fitstr, cix, true);
        set(fig, 'Position', [100, 100, 800, 1000]);
        fig_svfcn(fig, [dt '-' cellType '-' strjoin({num2str(cix)}, ',')]);
    end
end
close all

%%

figdir = 'figures/mt-hyperflow-overlays';
fig_fnfcn = @(tag) fullfile(figdir, [tag '-plain.png']);
fig_svfcn = @(fig, tag) hgexport(fig, fig_fnfcn(tag), ...
    hgexport('factorystyle'), 'Format', 'png');

fitdirf = @(mnkNm) ['data/evirepb-' mnkNm '/fits'];
fitstr = 'ASD';
cs = [2, 6, 1, 1, 6, 13, 14];
dts = {'20140304', '20140304', '20140305', '20140307', ...
    '20150324a', '20150324a', '20150324a'};
for ii = 1:numel(dts)
    if ii > 4
        isNancy = true;
        fitdir = fitdirf('nancy');
    else
        isNancy = false;
        fitdir = fitdirf('pat');
    end
    fig = plot.plotAllSaccadeKernelOverlays(dts{ii}, fitdir, isNancy, ...
        fitstr, cs(ii));
    axis off;
    fig_svfcn(fig, [dts{ii} '-MT-' num2str(cs(ii))]);
end
