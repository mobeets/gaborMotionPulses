function fitAllSTRFs(fitBehavior, fitCells, dts, mask, ...
    fitbasedir, figbasedir)
% fitSTRF(fitBehavior, fitCells, dts, mask, fitbasedir, figbasedir)
% 
% fits space-time receptive fields for each date in dts, using ASD and ML
%   if fitBehavior, fits monkey's psychophysical kernel
%   if fitCells, fits a cell's receptive field
%   if mask(1), fits using ASD
%   if mask(2), fits using ASD with grid search
%   if mask(3), fits using ML
%   if mask(4), fits using ASD with bilinear space and time weights
% 
% n.b. make sure to add path to github.com/mobeets/mASD
% 
    if nargin < 6
        figbasedir = 'figs-poiss';
    end
    if nargin < 5
        fitbasedir = 'fits-poiss';
    end
    if nargin < 4 || isempty(mask)
        mask = [true false false false]; % [ASD ASD_gs ML ASD_b]
    end
    if nargin < 3 || isempty(dts)
        dts = {'20130502', '20130514', '20130515', '20130517', ...
            '20130611', '20140213', '20140218', '20140226', '20140303', ...
            '20140304', '20140305', '20140306', '20140307', '20140310'};
    end
    if nargin < 2
        fitCells = true;
    end
    if nargin < 1
        fitBehavior = false;
    end
    
    if ~isempty(figbasedir) && ~exist(figbasedir, 'dir')
        mkdir(figbasedir);
    end
    if ~isempty(fitbasedir) && ~exist(fitbasedir, 'dir')
        mkdir(fitbasedir);
    end

    for ii = 1:numel(dts)
        dt = dts{ii};
        disp('----------------------');
        disp('----------------------');
        disp(dt);
        disp('----------------------');
        disp('----------------------');
        
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
        
        nfolds = 5;
        ifold = 1;
        devPct = 0.5;
        data = io.loadDataByDate(dt);
        % development/evaluation sets
        [~, evalinds] = reg.trainAndTest(data.X, data.R, devPct);
        [~, foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            llstr = 'poiss';
            scorestr = 'rsq';
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            
            lbs = [-3 -2 -1 -1]; ubs = [3 10 4 4]; ns = 7*ones(1,4);
            hyperOpts = struct('lbs', lbs, 'ubs', ubs, 'ns', ns, ...
                'isLog', true);
            Db = asd.sqdist.space(data.Xxy);
            
            % full ASD and ML
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr);
            BMAP = @(hyper) asd.fitHandle(hyper, Db, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);

            cell_inds = 1:size(data.Y_all, 2);
            ncells = numel(cell_inds);
            for nn = 1:ncells
                cell_ind = cell_inds(nn);
                lbl = [data.neurons{nn}.brainArea '-' num2str(cell_ind)];

                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                fits = fitSTRF(data, ML, MAP, BMAP, scoreFcn, ...
                    hyperOpts, figdir, lbl, ifold, mask, ...
                    foldinds, evalinds);

                if ~isempty(fitdir)
                    fits.isLinReg = true;
                    fits.foldinds = foldinds;
                    fits.evalinds = evalinds;
                    io.updateStruct(dat_fnfcn(lbl), fits);
                end

            end
        end
        %% run on decision
        if fitBehavior
            llstr = 'bern';
            scorestr = 'pseudoRsq';
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            
            lbs = [-3 -1 -1]; ubs = [3 4 4]; ns = 7*ones(1,3);
            hyperOpts = struct('lbs', lbs, 'ubs', ubs, 'ns', ns, ...
                'isLog', true);
            
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr);
            BMAP = @(hyper) asd.fitHandle(hyper, data.D, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);
            
            data.Y = data.R;

            lbl = 'decision';
            fits = fitSTRF(data, ML, MAP, BMAP, scoreFcn, hyperOpts, ...
                figdir, lbl, ifold, mask, foldinds, evalinds);
            if ~isempty(fitdir)
                fits.isLinReg = false;
                fits.foldinds = foldinds;
                fits.evalinds = evalinds;
                io.updateStruct(dat_fnfcn(lbl), fits);
            end
        end
        close all
    end
end
