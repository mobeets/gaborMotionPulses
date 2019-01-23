function [Y, X] = loadSpikeCounts(neurons, stim, binSizeSecs, kind)
% kind is one of ['motion', 'pulses']
    if nargin < 3
        binSizeSecs = 45/1000;
    end
    if nargin < 4
        kind = 'motion';
    end
    if strcmpi(kind, 'motion') % entire motion epoch
        tms = ([stim.timing.motionoff]+0.2) - [stim.timing.motionon];
        maxtm = max(tms);
        tbins = 0:binSizeSecs:maxtm;
        if maxtm ~= max(tbins)
            warning('Bin size will result in incomplete bins, which I ignore.');
        end
        nbins = numel(tbins);
    elseif strcmpi(kind, 'pulses') % motion epoch broken up by pulses
        nbins = 7;
    end
    
    X = nan(max(stim.trialnumber), nbins);
    Y = nan(max(stim.trialnumber), numel(neurons), nbins);
    for ii = 1:numel(neurons)
        neur = neurons{ii};
        if islogical(neur.trialIndex) && ~isequal(size(neur.trialIndex), ...
                size(neur.spikeCount)) || isempty(neur.trialIndex)
            continue;
        end
        sps = neur.spikeTimes;
        for jj = 1:numel(neur.trialIndex)
            ti = neur.trialIndex(jj);
            tmg = stim.timing(ti);
            t0 = tmg.plxstart;
            if strcmpi(kind, 'motion')
                t_starts = t0 + tmg.motionon + tbins;
                t_end = t0 + tmg.motionoff + 0.2;
                t_starts = t_starts(t_starts < t_end);
                t_ends = t_starts(2:end);
                t_starts = t_starts(1:end-1);                
            elseif strcmpi(kind, 'pulses')
                t_starts = t0 + tmg.pulses;
                t_ends = t_starts + 0.15 + 0.2;
            end
            for kk = 1:numel(t_starts)
                t0 = t_starts(kk);
                t1 = t_ends(kk);
                X(ti, kk) = kk;
                Y(ti, ii, kk) = sum(sps > t0 & sps <= t1);
            end
        end
    end
end
