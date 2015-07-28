function allSingleKernels(vs)
    ncells = numel(vs);
    aptrSz = 50;
    cntrSz = 300;    
    figure;
    
    % choose nrows, ncols for figure
    nr = 1; nc = ncells;
    if nc > 5
        nr = floor(sqrt(ncells));
        nc = ceil(ncells/nr);
    end
    
    % plot
    for ii = 1:ncells
        subplot(nr,nc,ii); hold on;
        v = vs(ii);
        if ~isfield(v, 'name')
            continue;
        end
        title(v.name, 'FontSize', 14);
        if isfield(v, 'wfSvd_1')
            plot.plotKernelSingle(v.Xxy, v.wfSvd_1(:,1), nan, aptrSz);
        else
            plot.plotKernelSingle(v.Xxy, v.wf(:,1), nan, aptrSz);
        end
        axis off;
        if isfield(v, 'rf_center0')
            scatter(v.rf_center0(:,1), v.rf_center0(:,2), cntrSz, 'r', ...
                'LineWidth', 2);
        end
    end
    
    width = 260;
    ny = width*nr;
    nx = 0.6*width*nc;
%     ny = 150; nx = 245;
    set(gcf, 'Position', [0 300 nx ny]);
end
