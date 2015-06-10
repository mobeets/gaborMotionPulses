function fitAllSTRFs2(runName, isNancy, fitType, dts)
% fitSTRF()
% 
% n.b. make sure to add path to github.com/mobeets/mASD
%
    if nargin < 4
        dts = {};
    end
    if isNancy
        mnkNm = 'nancy';
    else
        mnkNm = 'pat';
    end
    
    basedir = fullfile('data', [runName '-' mnkNm]);
    figbasedir = fullfile(basedir, 'figs');
    fitbasedir = fullfile(basedir, 'fits');
    
    fitCells = ~isempty(strfind(fitType, 'cell'));
    fitBehavior = ~isempty(strfind(fitType, 'behavior'));
    fitML = ~isempty(strfind(fitType, 'ML'));
    fitASD = ~isempty(strfind(fitType, 'ASD'));
    fitSpaceOnly = ~isempty(strfind(fitType, 'space-only'));
    if ~fitASD && ~fitML
        fitASD = true;
    end

    if isNancy && isempty(dts)
        dts = {'20150304a', '20150127' '20150304b', '20150305a', ...
            '20150305b', '20150305c', '20150306b', '20150306c', ...
            '20150309', '20150310', '20150312b', '20150313', ...
            '20150316c', '20150324a', '20150325', '20150326a', ...
            '20150331', '20150401', '20150402', '20150407a', ...
            '20150407b', '20150408a', '20150408b'};
    elseif isempty(dts)
        dts = {'20130502', '20130514', '20130515', '20130517', ...
            '20130611', '20140213', '20140218', '20140226', '20140303', ...
            '20140304', '20140305', '20140306', '20140307', '20140310'};
    end
    
    if ~isempty(basedir) && ~exist(basedir, 'dir')
        mkdir(basedir);
        mkdir(figbasedir);
        mkdir(fitbasedir);
    end
    
    warning off MATLAB:nearlySingularMatrix;
    for ii = 1:numel(dts)
        dt = dts{ii};
        disp('----------------------');
        disp(dt);
        
        if ~isempty(figbasedir)
            figdir = fullfile(figbasedir, dt);
            if ~exist(figdir, 'dir')
                mkdir(figdir);
            end
        else
            figdir = '';
        end
        if ~isempty(fitbasedir)
            fitdir = fullfile(fitbasedir, dt);
            if ~exist(fitdir, 'dir')
                mkdir(fitdir);
            end
        else
            fitdir = '';
        end
        dat_fnfcn = @(tag) fullfile(fitdir, [tag '.mat']);
                
        data = io.loadDataByDate(dt, isNancy);
        if fitSpaceOnly
            data.D = data.Ds;
            data.X = sum(data.Xf, 3);
            data.nt = 1;
        end
        % development/evaluation sets
        nfolds = 5;
        [~, foldinds] = tools.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            llstr = 'gauss';
            scorestr = 'rsq';

            cell_inds = 1:size(data.Y_all, 2);
            ncells = numel(cell_inds);
            for jj = 1:ncells
                cell_ind = cell_inds(jj);
                if data.neurons{cell_ind}.dPrime < 0.4
                    continue;
                end
                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                if sum(~isnan(data.Y)) == 0
                    continue;
                end
                label = [data.neurons{cell_ind}.brainArea '-' ...
                    num2str(cell_ind)];
                disp(label);

                if fitASD
                    obj = fitSTRF2(data, 'ASD', llstr, scorestr, ...
                        [label '-ASD'], foldinds);
                    plot.plotAndSaveKernel(obj, data, figdir);
                    fits.ASD = obj;
                end
                if fitML
                    obj = fitSTRF2(data, 'ML', llstr, scorestr, ...
                        [label '-ML'], foldinds);
                    plot.plotAndSaveKernel(obj, data, figdir);
                    fits.ML = obj;
                end
                if ~isempty(fitdir)
                    tools.updateStruct(dat_fnfcn(label), fits);
                end

            end
        end
        %% run on decision
        if fitBehavior
            llstr = 'bern';
            scorestr = 'pctCorrect';
            data.Y = data.R;
            label = 'decision';
            disp(label);
            
            if fitASD
                obj = fitSTRF2(data, 'ASD', llstr, scorestr, ...
                    [label '-ASD'], foldinds);
                plot.plotAndSaveKernel(obj, data, figdir);
                fits.ASD = obj;
            end
            if fitML
                obj = fitSTRF2(data, 'ML', llstr, scorestr, ...
                    [label '-ML'], foldinds);
                plot.plotAndSaveKernel(obj, data, figdir);
                fits.ML = obj;
            end
            if ~isempty(fitdir)
                tools.updateStruct(dat_fnfcn(label), fits);
            end
        end
        close all
    end
    warning on MATLAB:nearlySingularMatrix;
end
