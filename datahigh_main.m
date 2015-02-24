function D = datahigh_main(dt)
% dts = {'20130502', '20130514', '20130515', '20130517', '20130611', '20140213', '20140218', '20140226', '20140303', '20140304', '20140305', '20140306', '20140307', '20140310'};
    data = io.loadDataByDate(dt, '~/Desktop');
%     nsigfigs = nSigFigs(maxSpikeTime(data.neurons));
    nsigfigs = 3;
    D = spikesByTrial(data.stim, data.neurons, nsigfigs);
end

function maxTime = maxSpikeTime(neurons)
    maxTime = nan;
    for ii = 1:numel(neurons)
        maxTime = max(maxTime, max(neurons{ii}.spikeTimes));
    end 
end

function val = nSigFigs(maxTime)
    nr = maxTime - floor(maxTime);
    nrstr = num2str(nr,'%f');
    sigfgs = nrstr(nrstr~='0' & nrstr~='.');
    val = length(sigfgs);
end

function D = spikesByTrial(stim, neurons, nsigfigs)
    stimEventLength = 2.0;
    stimPreEvent = 0.2;
    inds = stim.goodtrial; % trials without broken fixation
    motionStartTimes = [stim.timing.motionon] - stimPreEvent + [stim.timing.plxstart];
    inds = inds & ~isnan(motionStartTimes)';
    ev1 = motionStartTimes(inds);
    ev2 = ev1 + stimEventLength; % const vector of motion length
    R = -(stim.targchosen(inds)-1) + 1;
    crct = stim.targcorrect == stim.targchosen;
    D = cell(numel(ev1), 3);
    for ii = 1:numel(ev1)
        D{ii,1} = spikesInWindow(neurons, ev1(ii), ev2(ii), nsigfigs)';
        D{ii,2} = 'traj';
        D{ii,3} = ['chc' num2str(R(ii))];% '-crct' num2str(crct(ii))]; % response is condition
    end
    D = cell2struct(D, {'data', 'type', 'condition'}, 2)';
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
