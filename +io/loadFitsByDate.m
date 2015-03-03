function fits = loadFitsByDate(dt, fitdir)
    fs = dir(fullfile(fitdir, dt));
    for ii = 1:numel(fs)
        f = fs(ii);
        if ~f.isdir && numel(f.name) > 3 && strcmp(f.name(end-2:end), 'mat')
            fname = strrep(f.name(1:end-4), '-', '_');
            fits.(fname) = load(fullfile(fitdir, dt, f.name));
        end
    end
end
