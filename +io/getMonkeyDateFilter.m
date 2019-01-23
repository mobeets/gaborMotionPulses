function ix = getMonkeyDateFilter(dts, mnks)
% gets dts mask for all monkey names present in mnks
    
    if isa(dts, 'cell') % if not already an array of numbers
        dts = cellfun(@(c) str2num(c(1:8)), dts);
    end
    ixP = dts < 20150101;
    ixN = ~ixP;
    
    % make mask
    ix = false(size(dts));
    if isempty(mnks) || isempty(mnks{1})
        % if empty, return true for all
        ix = ~ix;
        return;
    end
    ms = cellfun(@(m) m(1), mnks);
    if any(ismember(ms, 'P'))
        ix = ixP;
    end
    if any(ismember(ms, 'N'))
        ix = ix | ixN;
    end
end
