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
    [vals, data] = io.summaryByDate(dt, fitdir, 1);
    if nargin < 2 || any(isnan(cellind))
        cellinds = 1:numel(data.neurons);
    else
        cellinds = cellind;
    end
  
    figs = [];
    event = data.stim.targchosen;
    for ii = 1:numel(cellinds)
        [fig, name] = summaryBySingleCell(vals, data, event, cellinds(ii));
        if ~isempty(outdir)
            plot.saveFig(fig, name, outdir, figext);
        end
    end
end

function [fig, cellName] = summaryBySingleCell(vals, data, event, cellind)
    
    neuron = data.neurons{cellind};    
    [Z, bins] = io.psthByEvent(data.stim, neuron, event);
    
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
    
    subplot(2, 2, 1); hold on;
    lw = 1;
    for ii = 1:numel(Z)
        if ii == 1
            lbl = 'pref';
        else
            lbl = 'anti';
        end
        plot(bins(:,1)*1000, Z{ii}, '-', ...
            'Color', 'k', 'LineWidth', ii*lw, ...
            'DisplayName', lbl);
    end
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
    
    subplot(2, 2, 3); hold on;
    vmax = max(abs(wf(:)));
    plot.plotKernelSingle(data.Xxy, sgn*mu(:,1), vmax);
    title('spatial weights');
    xlabel(['sc=' fmt(sc) ', sep=' fmt(sep) ', r_{dec}=' fmt(r)]);
    
    subplot(2, 2, 4); hold on;
    bar(sgn*v(:,1), 'FaceColor', [0.3 0.8 0.3]);
    title('temporal weights');
    xlabel('pulse');

end
