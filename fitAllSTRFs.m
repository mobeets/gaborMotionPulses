function fitAllSTRFs(runName, ~, fitType, dts, brainArea)
% fitAllSTRFs(runName, ~, fitType, dts)
% 
% n.b. make sure to add path to mASD and cbrewer 
%
    if nargin < 4
        dts = {};
    end
    if nargin < 5
        brainArea = '';
    end
    if isempty(dts)
        stimdir = fullfile('data', 'stim');
        dts = io.getDates(stimdir);
%         stimdir = '..'; % assumes data in parent dir
%         dts = io.getDates(fullfile(stimdir, [mnkNm 'StimFiles']));
    end
    
    fitCells = ~isempty(strfind(fitType, 'cell'));
    fitBehavior = ~isempty(strfind(fitType, 'behavior'));
    fitML = ~isempty(strfind(fitType, 'ML'));
    fitASD = ~isempty(strfind(fitType, 'ASD'));
    fitFlat = ~isempty(strfind(fitType, 'Flat'));
    fitSpaceOnly = ~isempty(strfind(fitType, 'space-only'));
    if ~fitASD && ~fitML && ~fitFlat
        fitASD = true;
    end
    nfolds = 5;
    
    basedir = 'data';
    figbasedir = fullfile(basedir, 'figs', runName);
    fitbasedir = fullfile(basedir, 'fits', runName);
    if ~isempty(basedir) && ~exist(basedir, 'dir')
        mkdir(basedir);
        mkdir(figbasedir);
        mkdir(fitbasedir);
    end
    %%
    warning off MATLAB:nearlySingularMatrix;
    for ii = 1:numel(dts)
        dt = dts{ii};
        disp('----------------------');
        disp(dt);
        
        figdir = getOutputDir(figbasedir, dt);
        fitdir = getOutputDir(fitbasedir, dt);
                
        data = io.loadDataByDate(dt);
        if ~isfield(data, 'X')
            warning(['Error loading ' dt]);
            continue;
        end
        if size(data.X,1) < 5 % min # of trials
            continue;
        end
        if fitSpaceOnly
            data.D = data.Ds;
            data.X = sum(data.Xf, 3);
            data.nt = 1;
        end
        [~, foldinds] = tools.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            llstr = 'gauss';
            scorestr = 'rsq';

            cell_inds = 1:size(data.Y_all, 2);
            ncells = numel(cell_inds);
            for jj = 1:ncells
                cell_ind = cell_inds(jj);
%                 if data.neurons{cell_ind}.dPrime < 0.4
%                     continue;
%                 end
                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                if sum(~isnan(data.Y)) == 0
                    continue;
                end
                if ~isempty(brainArea) && ...
                        ~strcmpi(data.neurons{cell_ind}.brainArea, brainArea)
                    continue;
                end
                label = [data.neurons{cell_ind}.brainArea '-' ...
                    num2str(cell_ind)];
                disp(label);
                fitAndSaveSTRF(data, [fitASD fitML fitFlat], ...
                    {'ASD', 'ML', 'Flat'}, llstr, scorestr, label, ...
                    foldinds, figdir, fitdir);
            end
        end
        %% run on decision
        if fitBehavior
            llstr = 'bern';
            scorestr = 'mcc';
            data.Y = data.R;
            label = 'decision';
            disp(label);
            fitAndSaveSTRF(data, [fitASD fitML fitFlat], ...
                {'ASD', 'ML', 'Flat'}, llstr, scorestr, label, ...
                foldinds, figdir, fitdir);
        end
        close all
    end
    warning on MATLAB:nearlySingularMatrix;
end
%%
function outdir = getOutputDir(basedir, dt)
    if ~isempty(basedir)
        outdir = fullfile(basedir, dt);
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end
    else
        outdir = '';
    end
end
