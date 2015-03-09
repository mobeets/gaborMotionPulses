function dts = getDates(dirname)
    d = dir(dirname);
    folders = {d.name};
    dts = folders(cellfun(@(x) numel(x) > 2, folders));
end
