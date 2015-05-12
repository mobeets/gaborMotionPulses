function fig = kernelsByDate(dt, fitdir, isNancy)
    data = io.loadDataByDate(dt, isNancy);
    vals = io.summaryByDate(dt, data, fitdir, 1);
    dts = {vals.dt};
    isdt = strcmp(dts, dt);
    nms = {vals(isdt).name};
    typs = {vals(isdt).cellType};
    scs = {vals(isdt).score};
    seps = {vals(isdt).separability};
    yLblFcn = @(ii) [typs{ii} '-' nms{ii}(find(nms{ii}=='_')+1:end)];
    xLblFcn = @(ii) ['sc=' sprintf('%0.2f', scs{ii}) ...
        ', sep=' sprintf('%0.2f', seps{ii})];
    mus = horzcat(vals(isdt).mu0);
    vmx = repmat(max(abs(mus)), size(mus,1), 1);
    fig = plot.plotKernel(data.Xxy, mus./vmx, 1, nan, nan, 1.5, nan, ...
        xLblFcn, yLblFcn);
end
