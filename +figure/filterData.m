function vu = filterData(vu)
% filter fit data

    nMT = sum([vu.isMT]);

    % ignore sessions where monkey's pctCorrect < 75%
    vu = vu([vu.pctCorrect] >= 0.75);

    vuMT = vu([vu.isMT]);

    % ignore cells with dPrime < 0.4
%     vuMT = vuMT([vuMT.dPrime] >= 0.4);

    % ignore cells with ntrials < 100
    vuMT = vuMT([vuMT.ntrials] >= 100);

    % ignore deleted cells (just in case)
    badCells = {'20150407a_25', '20150407a_26', '20150407a_28', ...
        '20150407a_30', '20150407a_31', '20150407a_32', '20140304_14', ...
        '20140307_6', '20140307_8'};
    nms = arrayfun(@(x) [x.dt '_' num2str(x.id)], vuMT, 'uni', 0);
    ix = cellfun(@(c) sum(strcmp(c, badCells)) == 0, nms);
    assert(numel(ix) == numel(vuMT));
    vuMT = vuMT(ix);

    warning(['Removing ' num2str(nMT-numel(vuMT)) ' of ' ...
        num2str(nMT) ' MT cells.']);
    
    vu = [vu(~[vu.isCell]) vu([vu.isLIP]) vuMT];

end
