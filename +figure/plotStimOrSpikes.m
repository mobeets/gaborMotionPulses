function fig = plotStimOrSpikes(cells, cellName, startTrial, endTrial, kind)
% kind [str] is either 'stim' or 'spikes'

    cell = cells(strcmpi({cells.name}, cellName));
    trials = 1:numel(cell.ix);
    trials = trials(cell.ix);
    ix = (trials >= startTrial) & (trials <= endTrial);
    ts = trials(ix);
    X = cell.X(ix,:);
    y = cell.Y(ix);
    yh = cell.Yh(ix);
%     yh = cell.YhAR(ix);

    fig = plot.init(16);
    set(gca, 'LineWidth', 3);

    if strcmpi(kind, 'stim')
        % stimulus
        imagesc(ts, 1:size(X,2), X');
        colormap gray;
        set(gca, 'YTick', [1 size(X,2)]);
        ylim([1 size(X,2)]);
        ylabel('gabor index');
        
    else
        % spike counts (predicted vs. data)
        xx = [min(ts) ts max(ts)]; yy = [0 y' 0];
%         fill(xx(:), yy(:), 'k', 'FaceColor', 0.6*ones(3,1), 'Linestyle', 'none');
        plot(ts, y, 'LineWidth', 4, 'Color', 0.6*ones(3,1));
        plot(ts, yh, 'LineWidth', 4, 'Color', [0.2 0.45 0.2]);
        ylabel('spike count');
    end
    
    xlabel('trial');
    set(gca, 'XTick', [startTrial (startTrial+endTrial)/2 endTrial]);
    xlim([startTrial endTrial]);
    set(gca, 'TickDir', 'out');
    plot.setPrintSize(fig, struct('width', 9, 'height', 2.5));

end
