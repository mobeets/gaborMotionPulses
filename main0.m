%% load all

dts = {'20130502', '20130514', '20130515', '20130517', '20130611', '20140213', '20140218', '20140226', '20140303', '20140304', '20140305', '20140306', '20140307', '20140310'};
% dts = {'20130502'}
fitBehavior = false;
fitCells = true;

%% load

for ii = 1:numel(dts)
    dt = dts{ii};
    disp('----------------------');
    disp('----------------------');
    disp(dt);
    disp('----------------------');
    disp('----------------------');
    figdir = ['figs/' dt];
    datdir = ['fits/' dt];
    data = io.loadDataByDate(dt);
    % n.b. make sure to add path to github.com/mobeets/mASD

    %% params

    nfolds = 5;
    ifold = 1;
    [~,~,~,~,foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

    if ~exist(figdir, 'dir')
        mkdir(figdir);
    end
    if ~exist(datdir, 'dir')
        mkdir(datdir);
    end
    dat_fnfcn = @(tag) fullfile(datdir, [tag '.mat']);

    %% run on all cells
    if fitCells
        isLinReg = true;
        llstr = 'gauss';
        lbs = [-3, -2 -5 -5]; ubs = [3, 10 10 10]; ns = 5*ones(1,4);
        M = asd.linearASDStruct(data.D, llstr);
        mapFcn = M.mapFcn;
        scFcn = M.rsqFcn;
        mlFcn = @(~) ml.fitopts('gauss'); % no poisson for ML yet

    %     cell_inds = 2;
        cell_inds = 1:size(data.Y_all, 2);

        ncells = numel(cell_inds);
        for nn = 1:ncells
            cell_ind = cell_inds(nn);
            lbl = ['cell_' num2str(cell_ind)];

            data.Y = data.Y_all(:,cell_ind); % choose cell for analysis

            % ugly hack to drop nans in Y
            inds = ~isnan(data.Y);
            foldinds2 = foldinds(inds);
            data2 = data;
            data2.X = data2.X(inds,:,:);
            data2.Y = data2.Y(inds);
            fits = innermain(data2, foldinds2, mapFcn, mlFcn, scFcn, ...
                lbs, ubs, ns, figdir, lbl, ifold, [true false true]);

        %     fits = innermain(data, foldinds, mapFcn, mlFcn, scFcn, ...
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
        hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
        M = asd.logisticASDStruct(data.D);
        mapFcn = M.mapFcn;
        scFcn = M.pseudoRsqFcn;
        mlFcn = @(~) ml.fitopts('bern');
        data.Y = data.R;

        lbl = 'decision';
        fits = innermain(data, foldinds, mapFcn, mlFcn, scFcn, ...
            lbs, ubs, ns, figdir, lbl, ifold, [true false true]);
        if ~isempty(datdir)
            fits.isLinReg = isLinReg;
            io.updateStruct(dat_fnfcn(lbl), fits);
        end
    end
end
