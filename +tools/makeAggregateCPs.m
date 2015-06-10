function [cps, G] = makeAggregateCPs(vals, nm, grp)
    last_dt = '';
    Y = [];
    C = [];
    G = [];
    for ii = 1:numel(vals)
        if strcmp(vals(ii).type, 'decision')
            continue;
        end
        if ~strcmp(last_dt, vals(ii).dt)            
            d = io.loadDataByDate(vals(ii).dt);
            last_dt = vals(ii).dt;
        end        
        Y = [Y; vals(ii).(nm)];
        C = [C; vals(ii).C];
%         G = [G sum(d.X,2)];
%         G = [G; round(vals(ii).Ypos + vals(ii).Yneg)];
        G = [G; round(vals(ii).(grp))];
    end
    % conditional CP
    C = logical(C);
    [cps, G] = tools.AUC_singleByGroup(Y(C), Y(~C), ...
        unique(G), G(C), G(~C));
    
    figure; hold all;
    scatter(G, cps, 'filled');
    plot(xlim, [0.5 0.5], 'k--');
    xlabel(grp);
    ylabel(['CP(' nm ')']);
    title(['CP(' nm '), conditional on ' grp ...
        ', for all cells'' spikes combined']);
    ylim([0 1]);
end
