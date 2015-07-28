function objs = fitAllMethodsByName(name)    
    nfolds = 10;
    llstr = 'gauss';
    scorestr = 'rsq';
    figdir = '';
    fitdir = '';
    
    x = strsplit(name, '-');
    dt = x{1};
    x = strsplit(name, '_');
    cell_ind = str2num(x{2});
    
    d = io.loadDataByDate(dt);
    [~, foldinds] = tools.trainAndTestKFolds(d.X, d.R, nfolds);
    d.Y = d.Y_all(:, cell_ind);
    label = [d.neurons{cell_ind}.brainArea '-' num2str(cell_ind)];
    objs = fitAndSaveSTRF(d, [true true true], ...
        {'ASD', 'ML', 'Flat'}, llstr, scorestr, label, ...
        foldinds, figdir, fitdir);
    vs = struct();
    for ii = 1:numel(objs)
        v = objs{ii};
        vs(ii).Xxy = d.Xxy;
        vs(ii).name = [v.fitType ' (r^2=' sprintf('%0.2f', v.score) ')'];
        vs(ii).wf = reshape(v.w, v.shape(1), v.shape(2));
        [U, S, V] = svd(vs(ii).wf, 0);
        vs(ii).wfSvd_1 = U(:,1)*S(1)*V(:,1)';
    end
    plot.allSingleKernels(vs);
end
