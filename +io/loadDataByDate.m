function data = loadDataByDate(dt, isNancy, basedir, stimdir, ...
    spikesdir, ignoreFrozen, ignoreEarlyRepeats)
% 
    if nargin < 7
        ignoreEarlyRepeats = true;
    end
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
    Y = loadSpikeCounts(neurons, max(stim.trialnumber), stim.timing);
    frzinds = stim.frozentrials & stim.goodtrial;
    Yfrz = Y(frzinds,:);
    Rfrz = -(stim.targchosen(frzinds)-1) + 1;
    inds = stim.goodtrial;
     % each targ1 must be within 2 deg of median targ1, and same for targ2
%     targix = ignoreDistantTargets(stim, 2);
%     inds = inds & targix;
    if ignoreEarlyRepeats
        inds = inds & trialsWithoutRepeatStimStrengths(stim);
    end
    if ignoreFrozen
        inds = inds & ~stim.frozentrials;
    end
    Y = Y(inds,:);
    
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
    data.ixfrz = frzinds;
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

function Y = loadSpikeCounts(neurons, ny, stimtiming)
    Y = nan(ny, numel(neurons));
    Y2 = nan(ny, numel(neurons));
    for ii = 1:numel(neurons)
        cur = neurons{ii};
        if ~isequal(size(cur.trialIndex), size(cur.spikeCount))
            continue;
        end
%         sps = cur.spikeTimes;
%         for jj = 1:numel(cur.trialIndex)
%             ti = cur.trialIndex(jj);
%             tmg = stimtiming(ti);
%             t = tmg.plxstart;
%             t0 = t + tmg.motionon;
%             t1 = t + tmg.motionoff + 0.2;
%             Y2(ti, ii) = sum(sps >= t0 & sps <= t1);
%         end
        Y(cur.trialIndex, ii) = cur.spikeCount;
    end
%     Y = Y2;
end

function ix = trialsWithoutRepeatStimStrengths(stim)
    st = abs(sum(sum(stim.pulses, 3), 2));
    ix = true(size(st));
    if numel(st) < 100
        return;
    end
    st = st(1:100); % only consider first 100 trials
    us = unique(st(1:10));    
    if numel(us) < 3 % only two stim strengths used
        % look for repeats of these 1-2 strengths, ignore until past
        for ii = 1:numel(us)
            st(st==us(ii)) = inf;
        end
        ix0 = find(st ~= inf, 1);
        ix(1:ix0-1) = false;
    end    
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
