function data = loadDataByDate(dt, isNancy, basedir, stimdir, spikesdir, ignoreFrozen)    
    if nargin < 6
        ignoreFrozen = true;
    end    
    if nargin < 2
        isNancy = false;
    end
    if isNancy
        mnkNm = 'nancy';
    else
        mnkNm = 'pat';
    end
    if nargin < 5 || isempty(spikesdir)
        spikesdir = [mnkNm 'NeuronFiles'];
    end
    if nargin < 4 || isempty(stimdir)
        stimdir = [mnkNm 'StimFiles'];
    end
    if nargin < 3 || isempty(basedir)
        basedir = '/Volumes/LKCLAB/Users/Jay';
    end
    
    % load stimulus data
    stim = loadStim(dt, fullfile(basedir, stimdir));
    if isempty(fieldnames(stim))
        data = struct();
        return;
    end
    inds = stim.goodtrial;
    if ignoreFrozen
        inds = inds & ~stim.frozentrials;
    end
    X = stim.pulses;
    Xxy = stim.gaborXY;
    X = X(inds,:,:); % 1 (pref), -1 (anti-pref)
    [ny, nt, ns] = size(X);
    Xf = permute(X, [1 3 2]); % full stimulus: [ny x ns x nt]
    X = reshape(Xf, ny, ns*nt); % n.b. inverse is reshape(X, ny, ns, nt)
    D = asd.sqdist.spaceTime(Xxy, nt);
    Ds = asd.sqdist.space(Xxy);
    R = -(stim.targchosen(inds)-1) + 1; % 1->1 (pref), 2->0 (anti-pref)
    
    % load cell data
    neurons = loadNeurons(dt, fullfile(basedir, spikesdir));
    Y = loadSpikeCounts(neurons, max(stim.trialnumber));
    Y = Y(inds,:);
    
    % add all to output struct
    data.Xf = Xf;
    data.X = X;
    data.Y_all = Y;
    data.R = R;
    data.D = D;
    data.Ds = Ds;
    data.ndeltas = size(D, 3);
    data.Xxy = Xxy;
    data.ns = ns;
    data.nt = nt;
    data.stim = stim;
    data.neurons = neurons;
    data.dt = dt;
end

function stim = loadStim(dt, stimdir)
    stimfiles = io.findFile(stimdir, ['*' dt '_stim.mat'], true, true);
    if isempty(stimfiles)
        stim = struct();
        return;
    end
    stim = load(fullfile(stimdir, stimfiles{1}));
end

function neurons = loadNeurons(dt, spikesdir)
    spikefiles = io.findFile(spikesdir, ['*' dt '*.mat'], true, true);
    neurons = cell(numel(spikefiles), 1);
    for ii = 1:numel(spikefiles)
        neurons{ii} = load(fullfile(spikesdir, spikefiles{ii}));
    end
end

function Y = loadSpikeCounts(neurons, ny)
    Y = nan(ny, numel(neurons));
    for ii = 1:numel(neurons)
        cur = neurons{ii};
        if ~isequal(size(cur.trialIndex), size(cur.spikeCount))
            continue;
        end
        Y(cur.trialIndex, ii) = cur.spikeCount;
    end
end
