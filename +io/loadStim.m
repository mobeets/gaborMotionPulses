function stim = loadStim(dt, stimdir)
    stimfiles = io.findFile(stimdir, ['*' dt '_stim.mat'], true, true);
    if isempty(stimfiles)
        stim = struct();
        return;
    end
    stim = load(fullfile(stimdir, stimfiles{1}));
end
