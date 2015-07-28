function obs = fitAndSaveSTRF(data, fitMask, fitNames, llstr, scorestr, ...
    label, foldinds, figdir, fitdir, noPlot)

    dat_fnfcn = @(tag) fullfile(fitdir, [tag '.mat']);
    obs = cell(3,1);
    for ii = 1:numel(fitMask)
        if fitMask(ii)
            obj = fitSTRF(data, fitNames{ii}, llstr, scorestr, ...
                [label '-' fitNames{ii}], foldinds);
            plot.plotAndSaveKernel(obj, data, figdir, true, true, false);
            fits.(fitNames{ii}) = obj;
            obs{ii} = obj;
        end
    end
    if ~isempty(fitdir)
        tools.updateStruct(dat_fnfcn(label), fits);
    end
end
