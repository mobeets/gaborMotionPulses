function [objs, cps, ws] = simCPs(xopt, n)
    if nargin < 2
        n = 10;
    end
%     wL = linspace(-2, 0, n);
    wL = linspace(-0.1, 0, n);
%     wL = linspace(-0.8, -0.3, n);
    ws = [wL; ones(1, numel(wL))];
    objs = cell(numel(wL),1);
    cps = struct([]);
    wpos = cell(numel(wL),1);
    wneg = cell(numel(wL),1);
    for ii = 1:numel(wL)
        objs{ii} = sim.simCP(ws(:,ii), xopt, 233);
        cps = [cps objs{ii}.cp];
        wpos{ii} = objs{ii}.wpos;
        wneg{ii} = objs{ii}.wneg;
        % assert Ypos and Yneg correlation
    end
    
    wrs = sim.wposnegRatio(wpos, wneg);
%     figure; hist(wrs);
    figure, hold on; sim.plotCPs(cps, wrs, 'wr');
%     figure, hold on; sim.plotCPs(cps, [cps.cp_Y], 'CP(Y)');
    if xopt == 1
        titlestr = 'X perfectly correlated across space';
    elseif isnan(xopt)
        titlestr = 'X slightly correlated across space';
    else
        titlestr = 'X uncorrelated across space';
    end
    title(titlestr);

end
