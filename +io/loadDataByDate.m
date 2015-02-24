function data = loadDataByDate(dt, basedir, stimdir, spikesdir)    
    if nargin < 4 || isempty(spikesdir)
        spikesdir = 'SingleNeurons';
    end
    if nargin < 3 || isempty(stimdir)
        stimdir = 'stim';
    end
    if nargin < 2 || isempty(basedir)
        basedir = '~/Desktop';
    end
    
    % load stimulus data
    stim = loadStim(dt, fullfile(basedir, stimdir));    
    inds = stim.goodtrial;
    X = stim.pulses;
    Xxy = stim.gaborXY;
    X = X(inds,:,:); % 1 (pref), -1 (anti-pref)
    [ny, nt, ns] = size(X);
    X = permute(X, [1 3 2]);
    X = reshape(X, ny, nt*ns);
    D = asd.sqdist.spaceTime(Xxy, ns, nt);
    R = -(stim.targchosen(inds)-1) + 1; % 1->1 (pref), 2->0 (anti-pref)
    
    % load cell data
    neurons = loadNeurons(dt, fullfile(basedir, spikesdir));
    Y = loadSpikeCounts(neurons, max(stim.trialnumber));
    Y = Y(inds,:);
    
    % add all to output struct
    data.X = X;
    data.Y_all = Y;
    data.R = R;
    data.D = D;
    data.ndeltas = size(D, 3);
    data.Xxy = Xxy;
    data.ns = ns;
    data.nt = nt;
    data.stim = stim;
    data.neurons = neurons;
end

function stim = loadStim(dt, stimdir)
    stimfiles = io.findFile(stimdir, ['p' dt '_stim.mat'], true);
    stim = load(fullfile(stimdir, stimfiles{1}));
end

function neurons = loadNeurons(dt, spikesdir)
    spikefiles = io.findFile(spikesdir, dt);
    neurons = cell(numel(spikefiles), 1);
    for ii = 1:numel(spikefiles)
        neurons{ii} = load(fullfile(spikesdir, spikefiles{ii}));
    end
end

function Y = loadSpikeCounts(neurons, ny)
    Y = nan(ny, numel(neurons));
    for ii = 1:numel(neurons)
        cur = neurons{ii};
        Y(cur.trialIndex, ii) = cur.spikeCount;
    end
end
