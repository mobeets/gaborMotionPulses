function spatialWeightsByDate(vs, dt)

    tps = {'decision', 'MT', 'LIP'};
    vs0 = vs(strcmp({vs.dt}, dt));

    figure; hold on;
    nc = 5;
    nr = 4;
    c = 1;
    
    for jj = 1:numel(tps)
        vt0 = vs0(strcmp({vs0.type}, tps{jj}));
        if isempty(vt0)
            continue;
        end
        cmin = -max(arrayfun(@(x) max(abs(x.wfSvd_U(:,1))), vt0));
        for ii = 1:numel(vt0)
            subplot(nc, nr, c); hold on;
            vt0c = vt0(ii);
            w = vt0c.wfSvd_U(:,1);
            if sum(w) < 0
                w = -w;
            end
            plot.plotKernelSingle(vt0c.Xxy, w, nan, 50, cmin);
            
            val = vt0c.ecc;
            val2 = vt0c.dPrime;
            title([sprintf('%0.2f', val) ...
               ' = ' sprintf('%0.2f', val2)]);
%             title(vt0c.name);
            axis off;
            c = c + 1;
        end
    end
end
