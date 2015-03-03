function kernelsByDate(dt, fitdir)
    [vals, data] = io.loadSummariesByDate(dt, fitdir, 1);
    dts = {vals.dt};
    isdt = strcmp(dts, dt);
    nms = {vals(isdt).name};
    scs = {vals(isdt).score};
    seps = {vals(isdt).separability};
    lblFcn = @(ii) [nms{ii} ' (sc=' ...
        sprintf('%0.2f', scs{ii}) ', sep=' sprintf('%0.2f', seps{ii}) ')'];
    mus = horzcat(vals(isdt).mu0);
    vmx = repmat(max(abs(mus)), size(mus,1), 1);
    plot.plotKernel(data.Xxy, mus./vmx, nan, nan, nan, 1.5, nan, lblFcn);
end
