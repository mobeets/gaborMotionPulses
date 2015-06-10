isNancy = false;
fitdir = 'data/evirepb-pat/fits';
fitstr = 'ASD';

cellsPerFit = 3;

dts = io.getDates(fitdir);
for ii = 1:numel(dts)
    dt = dts{ii};
    d = io.loadDataByDate(dt, isNancy);
    inds = arrayfun(@(n) n.dPrime > 0.5 & ~isempty(n.hyperflow) & ...
        strcmp(n.brainArea, 'MT'), [d.neurons{:}]);
    ix = 1:numel(d.neurons);
    cellinds = ix(inds);
    nblks = floor(numel(cellinds)/cellsPerFit)+1;
    for jj = 1:nblks
        from = (jj-1)*cellsPerFit+1;
        to = ((jj-1)*cellsPerFit)+cellsPerFit;
        cix = cellinds(from:min(to,end));
        if isempty(cix)
            continue;
        end
        plot.plotAllSaccadeKernelOverlays(dt, fitdir, isNancy, fitstr, cix);
    end
end

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
