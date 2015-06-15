function fig = plotAllSaccadeKernelOverlays(dt, fitdir, isNancy, ...
    fitstr, cellinds, nSubplots)
% 
% Great MT Overlays:
%     20140304-2
%     20140304-6
%     20140305-1
%     20140307-1
% 
%     20150324a-6
%     20150324a-13
%     20150324a-14
% 
    if nargin < 6
        nSubplots = 1;
    end
    showPsth = true;
    ixf = @(jj, ii) (jj-1)*nSubplots+ii;
    
    d = io.loadDataByDate(dt, isNancy);
    event1 = d.stim.targchosen;
    event2 = sum(sum(d.stim.pulses, 3), 2) > 0; event2 = -(event2-2);
    fs = io.loadFitsByDate(dt, fitdir);
    if nargin < 5
        inds = arrayfun(@(n) n.dPrime > 0.5, [d.neurons{:}]);
        ix = 1:numel(d.neurons);
        cellinds = ix(inds);
%         cellinds = 1:numel(d.neurons);
    end
    fig = figure;
    for jj = 1:numel(cellinds)
        ii = cellinds(jj);
        n = d.neurons{ii};
        nm = [n.brainArea '_' num2str(ii)];
        if ~isfield(fs, nm)
            continue;
        end
        f = fs.(nm).(fitstr);
        subplot(numel(cellinds),nSubplots, ixf(jj,1)); hold on;
        plot.plotSaccadeKernelOverlay(d.stim, n, f{end}, true, true);
        
        if nSubplots > 1 && showPsth
            subplot(numel(cellinds), 2, ixf(jj,2)); hold on;            
            targPref = nanmax([n.targPref, 1]);
            if strcmp(n.brainArea, 'MT')
                event = event2;
                lbl = {'-mot', '+mot'};
            else
                event = event1;
                lbl = {'-chc', '+chc'};
            end
            plot.psthByEvent(d, n, event == targPref, ii, lbl);
            ylabel('spike rate by choice');
        elseif nSubplots > 1 && isfield(n, 'mtrfmap') && ~isempty(n.mtrfmap)
            subplot(numel(cellinds), ixf(jj,2)); hold all;
            thetas = [n.mtrfmap.thetas; n.mtrfmap.thetas(1)]/180*pi;
            rate   = [n.mtrfmap.rateMu; n.mtrfmap.rateMu(1)];
            polar(thetas, rate);
        end
    end
    
end
