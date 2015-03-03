function Z = psthByEvent(dt, eventName)
    if nargin < 2
        eventName = 'targchosen';
    end
    data = io.loadDataByDate(dt);
    event = data.stim.(eventName);
    
    figure; hold on;
    lw = 1;
    title(data.stim.exname);
    xlabel('time after motion onset (msec)');
    ylabel('average spike rate');
    
    clrs = lines(numel(data.neurons)+1);
    for ii = 1:numel(data.neurons)
        neuron = data.neurons{ii};
        subplot(numel(data.neurons), 1, ii); hold on;
        lbl = [neuron.brainArea '-' num2str(ii)];
        ylabel(lbl);
        
        [Z, bins] = io.psthByEvent(data.stim, neuron, event);
        for jj = 1:numel(Z)
            plot(bins(:,1)*1000, Z{jj}, repmat('-', 1, jj), ...
                'Color', clrs(ii,:), 'LineWidth', jj*lw, ...
                'DisplayName', lbl);
        end
    end
end
