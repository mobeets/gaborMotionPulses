function D = datahigh_main(dt)
    data = io.loadDataByDate(dt);
    D = spikesByTrial(data.stim, data.neurons, 3);
end

function D = spikesByTrial(stim, neurons, nsigfigs)
    stimEventLength = 2.0;
    stimPreEvent = 0.2;
    inds = stim.goodtrial; % trials without broken fixation
    motionStartTimes = [stim.timing.motionon] - stimPreEvent + [stim.timing.plxstart];
    inds = inds & ~isnan(motionStartTimes)';
    ev1 = motionStartTimes(inds);
    ev2 = ev1 + stimEventLength; % const vector of motion length
    dirprob = stim.dirprob(inds);
    chc = -(stim.targchosen(inds)-1) + 1;
    crct = stim.targcorrect == stim.targchosen;
    D = cell(numel(ev1), 5);
    for ii = 1:numel(ev1)
        D{ii,1} = spikesInWindow(neurons, ev1(ii), ev2(ii), nsigfigs)';
        D{ii,2} = 'traj';
        D{ii,3} = num2str(dirprob(ii));
        D{ii,4} = ['chc' num2str(chc(ii))];
        D{ii,5} = ['crct' num2str(crct(ii))];
    end
    D = cell2struct(D, ...
        {'data', 'type', 'condition', 'choice', 'correct'}, 2)';
end

function spikes = spikesInWindow(neurons, msStart, msEnd, nsigfigs)
    nbins = round((msEnd-msStart)*(10^nsigfigs))+1;
    nneurs = numel(neurons);
    spikes = zeros(nbins, nneurs);
    for ii = 1:nneurs
        sptms = neurons{ii}.spikeTimes;
        spiketimes = sptms(sptms >= msStart & sptms < msEnd) - msStart;
        inds = round(spiketimes*(10^nsigfigs))+1;
        sps = spikes(:,ii);
        sps(inds) = 1;
        spikes(:,ii) = sps;
    end
end
