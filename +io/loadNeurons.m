function neurons = loadNeurons(dt, spikesdir)
    spikefiles = io.findFile(spikesdir, ['*' dt '*.mat'], true, true);
    neurons = cell(numel(spikefiles), 1);
    for ii = 1:numel(spikefiles)
        neurons{ii} = load(fullfile(spikesdir, spikefiles{ii}));
    end
end