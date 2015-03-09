function fitAllSTRFs(fitBehavior, fitCells, dts)
% fitSTRF(fitBehavior, fitCells, dts)
% 
% fits space-time receptive fields for each date in dts, using ASD and ML
%   if fitBehavior, fits monkey's psychophysical kernel
%   if fitCells, fits a cell's receptive field
% 
% n.b. make sure to add path to github.com/mobeets/mASD
% 
    if nargin < 3
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
        figdir = ['figs/' dt];
        datdir = ['fits/' dt];                

        if ~exist(figdir, 'dir')
            mkdir(figdir);
        end
        if ~exist(datdir, 'dir')
            mkdir(datdir);
        end
        dat_fnfcn = @(tag) fullfile(datdir, [tag '.mat']);
        
        nfolds = 5;
        ifold = 1;
        data = io.loadDataByDate(dt);
        [~,~,~,~,foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            isLinReg = true;
            llstr = 'gauss';
            lbs = [-3, -2 -5 -5]; ubs = [3, 10 10 10]; ns = 5*ones(1,4);
            M = asd.linearASDStruct(data.D, llstr);
            mapFcn = M.mapFcn;
            scFcn = M.rsqFcn;
            mlFcn = @(~) ml.fitopts('gauss'); % no poisson for ML yet

            cell_inds = 1:size(data.Y_all, 2);
            ncells = numel(cell_inds);
            for nn = 1:ncells
                cell_ind = cell_inds(nn);
                lbl = [data.neurons{nn}.brainArea '-' num2str(cell_ind)];

                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis

                % ugly hack to drop nans in Y
                inds = ~isnan(data.Y);
                foldinds2 = foldinds(inds);
                data2 = data;
                data2.X = data2.X(inds,:,:);
                data2.Y = data2.Y(inds);
                fits = fitSTRF(data2, foldinds2, mapFcn, mlFcn, scFcn, ...
                    lbs, ubs, ns, figdir, lbl, ifold, [true false true]);

            %     fits = fitSTRF(data, foldinds, mapFcn, mlFcn, scFcn, ...
            %         lbs, ubs, ns, figdir, lbl, ifold, [true false true]);
                if ~isempty(datdir)
                    fits.isLinReg = isLinReg;
                    fits.foldinds = foldinds;
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

            lbl = 'decision';
            fits = fitSTRF(data, foldinds, mapFcn, mlFcn, scFcn, ...
                lbs, ubs, ns, figdir, lbl, ifold, [true false true]);
            if ~isempty(datdir)
                fits.isLinReg = isLinReg;
                io.updateStruct(dat_fnfcn(lbl), fits);
            end
        end
    end
end
