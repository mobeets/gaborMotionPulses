function dts = getDates(dirname)
    d = dir(dirname);
    folders = {d.name};
    dts = folders(cellfun(@(x) numel(x) > 2, folders));
    
    dts = dts(~strncmp(dts, '.', 1)); % ignore if starts with '.'
    
    % in case they aren't plain dates
    dts = strrep(dts, '_stim.mat', '');
    dts = strrep(dts, 'n', '');
    dts = strrep(dts, 'p', '');
end
