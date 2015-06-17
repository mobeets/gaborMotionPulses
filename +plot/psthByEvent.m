function psthByEvent(data, neuron, event, cellind, lbls)
    secsPerBin = 0.02;
    [Z, bins] = tools.getPsthByEvent(data.stim, neuron, event, 1.35, ...
        secsPerBin, 0.01);
    for ii = 1:numel(Z)
        if ii == 1
            lbl = lbls{ii}; % 'pref';
        else
            lbl = lbls{ii}; % 'anti';
        end
        plot(bins*1000, (1/secsPerBin)*Z{ii}, '-', ...
            'Color', 'k', 'LineWidth', ii, ...
            'DisplayName', lbl);
    end
    xlim([min(bins)*1000, max(bins)*1000]);
    legend('Location', 'NorthEastOutside');
    title([data.stim.exname '-' neuron.brainArea '-' num2str(cellind)]);
    xlabel('time after motion onset (msec)');
end
