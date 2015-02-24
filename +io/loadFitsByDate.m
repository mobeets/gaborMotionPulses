function fits = loadFitsByDate(dt, fitdir)
    fs = dir(fullfile(fitdir, dt));
    for ii = 1:numel(fs)
        f = fs(ii);
        if ~f.isdir && numel(f.name) > 3 && strcmp(f.name(end-2:end), 'mat')
            fname = f.name(1:end-4);
            fits.(fname) = load(fullfile(fitdir, dt, f.name));
        end
    end
end
