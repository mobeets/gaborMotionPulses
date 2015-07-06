function data = loadDataByDate2(dt, isNancy, basedir, stimdir, spikesdir, ignoreFrozen)    
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
%         basedir = '/Volumes/LKCLAB/Users/Jay';
        basedir = '/Users/jayhennig/Documents';
    end
    
    % load stimulus data
    stim = loadStim(dt, fullfile(basedir, stimdir));
    if isempty(fieldnames(stim))
        data = struct();
        disp(['ERROR: Could not find ' dt ' in ' fullfile(basedir, stimdir)]);
        return;
    end
    
    neurons = loadNeurons(dt, fullfile(basedir, spikesdir));
    [Y, sps] = loadSpikeCounts(neurons, max(stim.trialnumber), stim.timing);
    
    frzinds = stim.frozentrials & stim.goodtrial;
    Yfrz = Y(frzinds,:);
    Rfrz = -(stim.targchosen(frzinds)-1) + 1;
    inds = stim.goodtrial;
     % each targ1 must be within 2 deg of median targ1, and same for targ2
    targix = ignoreDistantTargets(stim, 2);
    inds = inds & targix;
    if ignoreFrozen
        inds = inds & ~stim.frozentrials;
    end
    Y = Y(inds,:);
    sps = sps(inds,:,:);
    data.sps = sps;
    
    X = stim.pulses;
    Xxy = stim.gaborXY;
    X = X(inds,:,:); % 1 (pref), -1 (anti-pref)
    [ny, nt, ns] = size(X);
    Xf = permute(X, [1 3 2]); % full stimulus: [ny x ns x nt]
    X = reshape(Xf, ny, ns*nt); % n.b. inverse is reshape(X, ny, ns, nt)
    D = asd.sqdist.spaceTime(Xxy, nt);
    Ds = asd.sqdist.space(Xxy);
    R = -(stim.targchosen(inds)-1) + 1; % 1->1 (pref), 2->0 (anti-pref)

    % add all to output struct
    data.ix = inds;
    data.Xf = Xf;
    data.X = X;
    data.Y_all = Y;
    data.Y_frz = Yfrz;
    data.R = R;
    data.R_frz = Rfrz;
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

function [Y, sps] = loadSpikeCounts(neurons, ny, stimtiming, prec)
    if nargin < 4
        prec = 20; % bins per second
    end
    Y = cell(ny, numel(neurons));
    for ii = 1:numel(neurons)
        n = neurons{ii};
        if ~isequal(size(n.trialIndex), size(n.spikeCount))
            continue;
        end
        sps = n.spikeTimes;
        for jj = 1:numel(n.trialIndex)
            ti = n.trialIndex(jj);
            tmg = stimtiming(ti);
            t = tmg.plxstart;
            t0 = t + tmg.motionon;
            t1 = t + tmg.motionoff;
            Y{ti,ii} = sps(sps >= t0 & sps <= t1) - t0;
        end
    end
    if nargout < 2
        return;
    end
    tmax = max(cell2mat(cellfun(@(x) max(x), Y(:), 'uni', 0)));
    bins = 0:(1/prec):ceil(tmax);
    nbins = numel(bins);
    sps = nan(ny, numel(neurons), nbins-1);
    for ii = 1:numel(neurons)
        for jj = 1:numel(n.trialIndex)
            ti = n.trialIndex(jj);
            Y0 = Y{ti,ii};
            for kk = 1:(nbins-1)
                sps(ti,ii,kk) = sum(Y0 >= bins(kk) & Y0 <= bins(kk+1));
            end
        end
    end
    Y = cellfun(@(x) numel(x), Y);
end

function ix = ignoreDistantTargets(stim, maxDist)
    ix1 = sqrt(sum((stim.targ1XY - repmat(median(stim.targ1XY), ...
        size(stim.targ1XY,1),1)).^2,2));
    ix2 = sqrt(sum((stim.targ2XY - repmat(median(stim.targ2XY), ...
        size(stim.targ1XY,1),1)).^2,2));
    ix = (ix1 <= maxDist) & (ix2 <= maxDist);
    
    if sum(~ix) > 0
        disp(['WARNING: removing ' num2str(sum(~ix)) ...
            ' trials due to varying target location.']);
    end
end
