function cellPairKernelsSingle(vu, n0, n1, dt)
    if isfloat(n0) % indexed by cellind and dt
        vut = vu(strcmp({vu.dt}, dt));
        v0 = vut([vut.cellind] == n0);
        v1 = vut([vut.cellind] == n1);
    else % indexed by name
        v0 = vu(strcmp({vu.name}, n0));
        v1 = vu(strcmp({vu.name}, n1));
    end
    aptrSz = 50;
    cntrSz = 300;
    
    figure;
    subplot(1,2,1); hold on;
    title(v0.name);
    plot.plotKernelSingle(v0.Xxy, v0.wfSvd_1(:,1), nan, aptrSz);
    axis off;  
    scatter(v0.rf_center0(:,1), v0.rf_center0(:,2), cntrSz, 'r', ...
        'LineWidth', 2);

    subplot(1,2,2); hold on;
    title(v0.name);
    plot.plotKernelSingle(v1.Xxy, v1.wfSvd_1(:,1), nan, aptrSz);
    axis off;
    scatter(v1.rf_center0(:,1), v1.rf_center0(:,2), cntrSz, 'r', ...
        'LineWidth', 2);
    set(gcf, 'Position', [0 300 200 150]);
    
end
