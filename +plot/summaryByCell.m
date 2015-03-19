function figs = summaryByCell(dt, cellind, fitdir, outdir, figext)
    if nargin < 5
        figext = 'png';
    end
    if nargin < 4
        outdir = '';
    end
    if nargin < 3
        fitdir = 'fits';
    end
    data = io.loadDataByDate(dt);
    vals = io.summaryByDate(dt, data, fitdir, 1);
    if nargin < 2 || any(isnan(cellind))
        cellinds = 1:numel(data.neurons);
    else
        cellinds = cellind;
    end
  
    figs = [];
    event1 = data.stim.targchosen;
    event2 = sum(sum(data.stim.pulses, 3), 2) > 0; % 0->neg, 1->pos
    for ii = 1:numel(cellinds)
        [fig, name] = summaryBySingleCell(dt, vals, data, ...
            event1, event2, cellinds(ii));
        figs = [figs fig];
        if ~isempty(outdir)
            plot.saveFig(fig, name, outdir, figext);
        end
    end
end

function [fig, cellName] = summaryBySingleCell(dt, vals, data, event1, ...
    event2, cellind)
    
    neuron = data.neurons{cellind};
    targPref = nanmax([neuron.targPref, 1]);
    
    cellName = [neuron.brainArea '_' num2str(cellind)];
    ind = strcmp({vals.name}, cellName);
    wf = vals(ind).mu;
    [u,s,v] = svd(wf);
    sgn = sign(sum(v(:,1)));
    mu = u(:,1)*s(1)*v(:,1)';
    sc = vals(ind).score;
    sep = vals(ind).separability;
    r = vals(ind).decisionCorrelation;
    fmt = @(val) sprintf('%0.2f', val);
    
    fig = figure;
    set(gcf,'color','w');
    pos1 = get(fig, 'Position');
    new_pos1 = pos1.*[1 1 1 1.5];
    set(fig, 'Position', new_pos1);
    
    subplot(3, 2, 1); hold on;
    [Z1, b1] = io.getPsthByEvent(data.stim, neuron, ...
        event1 == targPref);
    plotPsth(Z1, b1, data, neuron, cellind, {'-chc', '+chc'});
    ylabel('spike rate by choice');
    
    subplot(3, 2, 2); hold on;
    [Z2, b2] = io.getPsthByEvent(data.stim, neuron, ...
        -(event2-2) == targPref);
    plotPsth(Z2, b2, data, neuron, cellind, {'-mot', '+mot'});
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

function plotPsth(Z, bins, data, neuron, cellind, lbls)
    lw = 1;
    for ii = 1:numel(Z)
        if ii == 1
            lbl = lbls{ii}; % 'pref';
        else
            lbl = lbls{ii}; % 'anti';
        end
        plot(bins*1000, Z{ii}, '-', ...
            'Color', 'k', 'LineWidth', ii*lw, ...
            'DisplayName', lbl);
    end
    xlim([min(bins)*1000, max(bins)*1000]);
    legend('Location', 'NorthEastOutside');
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
end
