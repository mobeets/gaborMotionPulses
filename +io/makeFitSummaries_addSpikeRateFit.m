function makeFitSummaries_addSpikeRateFit(vs)
    lastdt = nan;
    for ii = 1:numel(vs)
        v = vs(ii);    
        if isnan(v.cellind)
            vs(ii).spikeRateSlope = nan;
            vs(ii).spikeRateRsq = nan;
            continue;
        end
        if ~strcmp(v.dt, lastdt)
            close all;
            lastdt = v.dt;
            [ps, nms] = plot.spikeRateByTrial(lastdt); 
        end
        nm = [v.dt '-' num2str(v.cellind)];
        p = ps(strcmp(nms, nm),:);
        vs(ii).spikeRateSlope = p(1);
        vs(ii).spikeRateRsq = p(3);
    end
end
