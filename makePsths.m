%% load

dt = 'n20150324a';
basedir = 'data';
spikesdir = 'neurons';
stimdir = 'stim';
% binSizeSecs = 45/2000;
nPulses = 7;
nBinsPerPulse = 7;
binSizeSecs = 0.15/nBinsPerPulse; % to evenly divide pulse length

stim = io.loadStim(dt, fullfile(basedir, stimdir));
neurons = io.loadNeurons(dt, fullfile(basedir, spikesdir));
[Y, Bins] = io.loadSpikeCounts(neurons, stim, binSizeSecs);
disp('Done');

%% find psths and plot

ix = stim.goodtrial & stim.correct;
Xdirprob = stim.dirprob(ix,:);
Pulses = stim.pulses(ix,:,:);
Yall = Y(ix,:,:);
Yall = Yall(:,:,1:55);
warning('Chopping off last few time steps.');

% split up by pulse inds and flatten
pulseInds = (1:nBinsPerPulse:nPulses*nBinsPerPulse) + 4;
pulseInds = pulseInds(1);
pulseInds = 1;
% pulseInds = pulseInds(2:(end-1)); % only keep 2nd-6th pulses

pulseLocationsToAvg = 1:size(Pulses,3); % which pulse locations to avg over
% pulseLocationsToAvg = 1:7;
% pulseLocationsToAvg = (19-7+1):19;

Xpulse = [];
Ypulse = [];
for ii = 1:numel(pulseInds)
    t1 = pulseInds(ii);
    t2 = pulseInds(ii) + nBinsPerPulse + 4;
    yc = Yall(:,:,t1:t2);
    xc = sum(Pulses(:,ii,pulseLocationsToAvg),3);
%     xc = unique(round(xc / 4)); % discretize x
    Xpulse = [Xpulse; xc];
    Ypulse = [Ypulse; yc];
end

%% make psths

% Xc = Xdirprob; Yc = Yall;
Xc = Xpulse; Yc = Ypulse;
xs = unique(Xc);
xs = xs(3:(end-2));

psths = nan(numel(xs), size(Yc,2), size(Yc,3));
for ii = 1:size(xs,1)
    psths(ii,:,:) = nanmean(Yc(Xc == xs(ii),:,:),1);
end

pulseToHighlight = nan;
if ~isnan(pulseToHighlight)
    pulseInds = (1:nBinsPerPulse:nPulses*nBinsPerPulse) + 4;
    highlightInds = pulseInds(pulseToHighlight):pulseInds(pulseToHighlight+1);
else
    highlightInds = [];
end

clrs = cbrewer('div', 'RdBu', numel(xs));
plot.init;
for ii = 1:size(psths,2)
    subplot(4,4,ii); hold on;    
    for jj = 1:size(psths,1)
        yc = squeeze(psths(jj,ii,:));
        if isempty(highlightInds)
            plot(yc, '-', 'Color', clrs(jj,:));
        else
            plot(yc, '-', 'Color', 0.8*ones(3,1));
            plot(highlightInds, yc(highlightInds), '-', ...
                'LineWidth', 2, 'Color', clrs(jj,:));
        end        
    end
    xlim([0 size(psths,3)]);
end

%% apply pca and plot

% Xc = Xdirprob; Yc = Yall;
Xc = Xpulse; Yc = Ypulse;

Yflat = nan(size(Yc,1)*size(Yc,3), size(Yc,2));
for ii = 1:size(Yc,2)
    yc = Yc(:,ii,:);
    Yflat(:,ii) = yc(:);
end

X2 = repmat(Xc,1,size(Yc,3))';
tms = repmat(1:size(Yc,3),1,size(Xc,1));
Xs = [X2(:) tms'];

[coeff, score, ~, ~, vexp, mu] = pca(Yflat);

% show dimensionality
plot.init;
plot(cumsum(vexp), '.-');
xlabel('# PCs');
ylabel('Cumulative % Variance Explained');
ylim([0 100]);

% plot psths of PC's
plot.init;
for ii = 1:numel(xs)
    pstc = squeeze(psths(ii,:,:))';
    yc = bsxfun(@minus, pstc, mu)*coeff;
    for jj = 1:size(yc,2)
        subplot(4,4,jj); hold on;
        plot(binSizeSecs*(1:size(yc,1)), yc(:,jj), '-', 'Color', clrs(ii,:));
        xlim([0 binSizeSecs*size(yc,1)]);
    end
end

%% psth movie

% movie of timecourse of psths in 3D
nreps = 5;
vmx = max(abs(psths(:)));
plot.init;
for jj = 1:nreps
    for t = 1:size(psths,3)
        for ii = 1:numel(xs)
            pstc = squeeze(psths(ii,:,:))';
            yc = bsxfun(@minus, pstc, mu)*coeff;
            plot3(yc(1:t,1), yc(1:t,2), yc(1:t,3), '.-', 'Color', clrs(ii,:));
        end
        xlim([-vmx vmx]); ylim(xlim); zlim(xlim);
%         axis equal;
%         view(-40, 70);
        pause(0.2);
        title(['t = ' num2str(t*binSizeSecs)]);
    end
    if jj < nreps
        cla;
    end
end
