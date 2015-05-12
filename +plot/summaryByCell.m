function figs = summaryByCell(dt, cellind, isNancy, fitdir, outdir, figext)
    if nargin < 6
        figext = 'png';
    end
    if nargin < 5
        outdir = '';
    end
    if nargin < 4
        fitdir = 'fits';
    end
    data = io.loadDataByDate(dt, isNancy);
    vs = io.makeFitSummaries(fitdir, isNancy, 'ASD', {dt});
    if nargin < 2 || any(isnan(cellind))
        cellinds = 1:numel(data.neurons);
    else
        cellinds = cellind;
    end
  
    figs = [];
    event1 = data.stim.targchosen;
    event2 = sum(sum(data.stim.pulses, 3), 2) > 0; % 0->neg, 1->pos
    for ii = 1:numel(cellinds)
        [fig, name] = summaryBySingleCell(vs, data, ...
            event1, event2, cellinds(ii));
        figs = [figs fig];
        if ~isempty(outdir)
            plot.saveFig(fig, name, outdir, figext);
        end
    end
end

function [fig, cellName] = summaryBySingleCell(vs, data, event1, ...
    event2, cellind)
    
    neuron = data.neurons{cellind};
    targPref = nanmax([neuron.targPref, 1]);
    
    cellName = [neuron.brainArea '_' num2str(cellind)];
    ind = strcmp({vs.type}, neuron.brainArea) & ...
        ([vs.cellind] == cellind);
    wf = vs(ind).mu;
    [u,s,v] = svd(wf);
    sgn = sign(sum(v(:,1)));
    mu = u(:,1)*s(1)*v(:,1)';
    sc = vs(ind).score;
    sep = vs(ind).separability;
    r = vs(ind).decisionCorrelation;
    fmt = @(val) sprintf('%0.2f', val);
    
    fig = figure;
    set(gcf,'color','w');
    pos1 = get(fig, 'Position');
    new_pos1 = pos1.*[1 1 1 1.5];
    set(fig, 'Position', new_pos1);
    
    subplot(3, 2, 1); hold on;
    plot.psthByEvent(data, neuron, event1 == targPref, ...
        cellind, {'-chc', '+chc'});
    ylabel('spike rate by choice');
    
    subplot(3, 2, 2); hold on;
    plot.psthByEvent(data, neuron, -(event2-2) == targPref, ...
        cellind, {'-mot', '+mot'});
    ylabel('spike rate by motion dir');
    
    subplot(3, 2, 3); hold on;
    plot.plotKernelSingle(data.Xxy, sgn*mu(:,1), max(abs(mu(:))));
    title('spatial weights');
    xlabel(['sc=' fmt(sc) ', sep=' fmt(sep) ', r_{dec}=' fmt(r)]);
    
    subplot(3, 2, 4); hold on;
    bar(sgn*v(:,1), 'FaceColor', [0.3 0.8 0.3]);
    title('temporal weights');
    xlabel('pulse');
    
    subplot(3, 2, 5); hold on;
    [cp, xs] = io.getCP(data.stim, neuron, event1 == targPref, ...
        0.0, 1.5, 0.2, 0.1);
    plot(xs, cp);
    ylim([0.0 1.0]);
    title('CP');
    xlabel('time after motion onset (sec)');

end
