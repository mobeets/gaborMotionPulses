function [X, Y, Xxy, R] = makeDesign(infile, outfile, neuronID, stimEventLength)
% creates a design matrix of pulses (X) and spike counts (Y)
% infile must be .mat containing stim and spikes data
% if neuronID is int, loads that neuron's spike counts only;
% otherwise, loads all neurons' spike counts
%
% Outputs:
%   X [n x nw] - stimulus for nw gabors
%   Y [n x nc] - spike counts for nc cells
%   Xxy [nw x 2] - gabor xy-coordinates
%   R [n x 1] - monkey choice
% 
% requirements: https://github.com/jcbyts/pdstools
% 
    import pdsa.countSpikes
    if nargin < 1
        infile = 'testData';
    end
    if nargin < 2
        outfile = '';
    end
    if nargin < 3
        neuronID = 0; % load all neurons
    end
    if nargin < 4
        stimEventLength = 1.5;
    end
    load(infile, 'stim', 'spikes');
    
    %% load stimulus and timings
    inds = stim.goodtrial; % trials without broken fixation
    motionStartTimes = stim.timing.motionon(:,1) + stim.timing.plxstart;
    ev1 = motionStartTimes(inds);
    ev2 = ev1 + stimEventLength; % const vector of motion length
    X = stim.pulses(inds,:,:);
    R = -(stim.targchosen(inds)-1) + 1; % target choice; 2 -> 0, 1 -> 1
    Xxy = stim.gaborXY;

    %% load spike counts
    if neuronID > 0
        neuronIDs = neuronID;
    elseif isfield(spikes, 'goodinds')
        neuronIDs = spikes.goodinds;
    else
        neuronIDs = unique(spikes.id);
    end
    Y = nan(sum(inds), numel(neuronIDs));
    for ii = 1:numel(neuronIDs)
        sps = spikes.time(spikes.id==neuronIDs(ii));
        Y(:,ii) = countSpikes(sps, ev1, ev2);
    end

    %% save design mat
    if ~isempty(outfile)
        save(outfile, 'X', 'Y', 'Xxy', 'R');
    end
end
