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
        figbasedir = 'figs';
    end
    if nargin < 5
        fitbasedir = 'fits';
    end
    if nargin < 4 || isempty(mask)
        mask = [true false true true]; % [ASD ASD_gs ML ASD_b]
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
        fitBehavior = true;
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
        data = io.loadDataByDate(dt);
        [~,~,~,~,foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            isLinReg = true;
            llstr = 'gauss';
            lbs = [-3, -2 -5 -5]; ubs = [3, 10 10 10]; ns = 5*ones(1,4);
            
            % full ASD and ML
            M = asd.linearASDStruct(data.D, llstr);
            mapFcn = M.mapFcn;
            scFcn = M.rsqFcn;
            mlFcn = @(~) ml.fitopts('gauss'); % no poisson for ML yet
            
            % bilinear ASD
            opts = struct(); opts.shape = {data.ns, data.nt};
            M2 = asd.linearASDStruct(data.Ds, llstr, 'bilinear', opts);
            bmapFcn = M2.mapFcn;

            cell_inds = 1:size(data.Y_all, 2);
            ncells = numel(cell_inds);
            for nn = 1:ncells
                cell_ind = cell_inds(nn);
                lbl = [data.neurons{nn}.brainArea '-' num2str(cell_ind)];

                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                [cur_data, cur_foldinds] = dropTrialsIfNaN(data, ...
                    foldinds, ~isnan(data.Y));
                fits = fitSTRF(cur_data, cur_foldinds, mapFcn, mlFcn, ...
                    bmapFcn, scFcn, lbs, ubs, ns, figdir, lbl, ifold, ...
                    mask);

                if ~isempty(fitdir)
                    fits.isLinReg = isLinReg;
                    fits.foldinds = cur_foldinds;
                    io.updateStruct(dat_fnfcn(lbl), fits);
                end

            end
        end
        %% run on decision
        if fitBehavior
            isLinReg = false;
            lbs = [-3, -5 -5]; ubs = [3, 10 10]; ns = 5*ones(1,3);
            M = asd.logisticASDStruct(data.D);
            mapFcn = M.mapFcn;
            scFcn = M.pseudoRsqFcn;
            mlFcn = @(~) ml.fitopts('bern');
            data.Y = data.R;
            
            % bilinear ASD
            opts = struct(); opts.shape = {data.ns, data.nt};
            M2 = asd.logisticASDStruct(data.Ds, 'bilinear', opts);
            bmapFcn = M2.mapFcn;

            lbl = 'decision';
            fits = fitSTRF(data, foldinds, mapFcn, mlFcn, bmapFcn, ...
                scFcn, lbs, ubs, ns, figdir, lbl, ifold, mask);
            if ~isempty(fitdir)
                fits.isLinReg = isLinReg;
                io.updateStruct(dat_fnfcn(lbl), fits);
            end
        end
    end
end

function [cur_data, cur_foldinds] = dropTrialsIfNaN(data, inds, foldinds)
    cur_foldinds = foldinds(inds);
    cur_data = data;
    cur_data.X = cur_data.X(inds,:);
    cur_data.Xf = cur_data.Xf(inds,:,:);
    cur_data.Y = cur_data.Y(inds);
    cur_data.R = cur_data.R(inds);
end

