function fig = deltaDecodingAcc(pairs)

    % compute score change from shuffle
    scsDeltaMean = [pairs.scoreGainWithCorrs];
    scsDeltaSe = [pairs.scoreGainWithCorrs_se];
    scsDelta = [scsDeltaMean-scsDeltaSe; scsDeltaMean; scsDeltaMean+scsDeltaSe]';

    fig = plot.init(18);
    set(gca, 'LineWidth', 2);
    [~,ix] = sort(scsDelta(:,2));
    plot([1 size(scsDelta,1)], [0 0], '-', 'LineWidth', 2, ...
        'Color', 0.8*ones(3,1));
    plot(100*scsDelta(ix,:), '-', 'LineWidth', 2, ...
        'Color', [0.4 0.4 0.4]);
    xlabel('pair index (sorted)');
    ylabel('\Delta decoding accuracy');
    % ylabel({'\Delta decoding accuracy', ...
    %     '\leftarrow Correlations hurt      Correlations help \rightarrow'});
    ytcks = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(n) [num2str(n) '%'], ytcks, 'uni', 0));
%     set(gca, 'TickLength', [0 0]);
    plot.setPrintSize(gcf, struct('width', 4, 'height', 3));

end