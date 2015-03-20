function psthByEvent(data, neuron, event, cellind, lbls)
    [Z, bins] = io.getPsthByEvent(data.stim, neuron, event);
    for ii = 1:numel(Z)
        if ii == 1
            lbl = lbls{ii}; % 'pref';
        else
            lbl = lbls{ii}; % 'anti';
        end
        plot(bins*1000, Z{ii}, '-', ...
            'Color', 'k', 'LineWidth', ii, ...
            'DisplayName', lbl);
    end
    xlim([min(bins)*1000, max(bins)*1000]);
    legend('Location', 'NorthEastOutside');
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
end
