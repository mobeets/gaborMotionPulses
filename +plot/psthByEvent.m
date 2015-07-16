function psthByEvent(data, neuron, event, cellind, lbls, secsPerBin)
    if nargin < 6
        secsPerBin = 0.035;
    end
    [Z, bins] = tools.getPsthByEvent(data.stim, neuron, event, 1.35, ...
        secsPerBin, 0.01);
    for ii = 1:numel(Z)
        plot(bins*1000, (1/secsPerBin)*Z{ii}, '-', ...
            'Color', 'k', 'LineWidth', ii, ...
            'DisplayName', lbls{ii});
    end
    xlim([min(bins)*1000, max(bins)*1000]);
    legend('Location', 'NorthEastOutside');
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
end
