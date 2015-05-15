function fitAllSTRFs(runName, isNancy, fullTemporalSmoothing)
% fitSTRF()
% 
% n.b. make sure to add path to github.com/mobeets/mASD
%
    if isNancy
        mnkNm = 'nancy';
    else
        mnkNm = 'pat';
    end
    figbasedir = ['data/figs-' mnkNm '-' runName];
    fitbasedir = ['data/fits-' mnkNm '-' runName];
    opts.fullTemporalSmoothing = fullTemporalSmoothing;
        
%     figbasedir = 'data/figs-evi-nancy';
%     fitbasedir = 'data/fits-evi-nancy';
    fitMask = [false false false false true]; % [ASD_g ASD_gs ML ASD_b ASD]
    fitCells = true;
    fitBehavior = false;
    nfolds = 5;
    devPct = 1.0; % pct of data used to generate fit, after choosing hyper
    lbs = [-3 -2 -1 -1]; ubs = [3 10 6 6]; ns = 10*ones(1,4);
%     lbs = [0 0 1 1]; ubs = [0 0 1 1]; ns = 1*ones(1,4);    
%     lbs(end) = 5e3; ubs(end) = 5e3; ns(end) = 1; % fix delta_t
    
    if isNancy
        dts = {'20150304a', '20150127' '20150304b', '20150305a', ...
            '20150305b', '20150305c', '20150306b', '20150306c', ...
            '20150309', '20150310', '20150312b', '20150313', ...
            '20150316c', '20150324a', '20150325', '20150326a', ...
            '20150331', '20150401', '20150402', '20150407a', ...
            '20150407b', '20150408a', '20150408b'};
    else
        dts = {'20130502', '20130514', '20130515', '20130517', ...
            '20130611', '20140213', '20140218', '20140226', '20140303', ...
            '20140304', '20140305', '20140306', '20140307', '20140310'};
    end
    dts = {'20150304a'};
    
    if ~isempty(figbasedir) && ~exist(figbasedir, 'dir')
        mkdir(figbasedir);
    end
    if ~isempty(fitbasedir) && ~exist(fitbasedir, 'dir')
        mkdir(fitbasedir);
    end
    
    warning off MATLAB:nearlySingularMatrix;
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
                
        data = io.loadDataByDate(dt, isNancy);
        % development/evaluation sets
        [~, evalinds] = reg.trainAndTest(data.X, data.R, devPct);
        [~, foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            llstr = 'gauss';
            fitstr = 'evi';
            
            if strcmpi(llstr, 'gauss')                
                lbsc = lbs; ubsc = ubs; nsc = ns;                
                scorestr = 'rsq';
            else % remove ssq hyper
                scorestr = 'pseudoRsq';
                ix = [1 3 4]; lbsc = lbs(ix); ubsc = ubs(ix); nsc = ns(ix);
            end
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            hyperOpts = struct('lbs', lbsc, 'ubs', ubsc, 'ns', nsc, ...
                'isLog', true);
            
            % full ASD and ML
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr, fitstr, opts);
            Db = asd.sqdist.space(data.Xxy);
            BMAP = @(hyper) asd.fitHandle(hyper, Db, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);

            cell_inds = 1:size(data.Y_all, 2);
%             cell_inds = [2];
            ncells = numel(cell_inds);
            for nn = 1:ncells
                cell_ind = cell_inds(nn);
                lbl = [data.neurons{cell_ind}.brainArea '-' num2str(cell_ind)];

                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                if sum(~isnan(data.Y)) == 0
                    continue;
                end
                fits = fitSTRF(data, ML, MAP, BMAP, scoreFcn, ...
                    hyperOpts, figdir, lbl, fitMask, ...
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
            
            if strcmpi(llstr, 'gauss') || strcmp(llstr, 'evi')
                lbsc = lbs; ubsc = ubs; nsc = ns;
                scorestr = 'rsq';
            else % remove ssq hyper
                scorestr = 'pseudoRsq';
                ix = [1 3 4]; lbsc = lbs(ix); ubsc = ubs(ix); nsc = ns(ix);
            end
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            hyperOpts = struct('lbs', lbsc, 'ubs', ubsc, 'ns', nsc, ...
                'isLog', true); 
            
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr);
            BMAP = @(hyper) asd.fitHandle(hyper, data.D, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);
            data.Y = data.R;

            lbl = 'decision';
            fits = fitSTRF(data, ML, MAP, BMAP, scoreFcn, hyperOpts, ...
                figdir, lbl, fitMask, foldinds, evalinds);
            if ~isempty(fitdir)
                fits.isLinReg = false;
                fits.foldinds = foldinds;
                fits.evalinds = evalinds;
                io.updateStruct(dat_fnfcn(lbl), fits);
            end
        end
        close all
    end
    warning on MATLAB:nearlySingularMatrix;
end
