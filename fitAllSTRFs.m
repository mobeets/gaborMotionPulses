function fitAllSTRFs(runName, isNancy, fitType, dts)
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
    
    fitCells = strcmpi(fitType, 'cells');
    fitBehavior = strcmpi(fitType, 'behavior');
    fitMask = [false false false false false]; % [ASD_g ASD_gs ML ASD_b ASD]
    fitMask(1) = fitBehavior;
    fitMask(5) = fitCells;
    
    fitSpaceOnly = false;

    nfolds = 5;
    devPct = 1.0; % pct of data used to generate fit, after choosing hyper
    lbs = [-10 -5 -1 -1]; ubs = [10 6 6 6]; ns = 7*ones(1,4);
    if fitSpaceOnly
        lbs = lbs(1:end-1); ubs = ubs(1:end-1); ns = ns(1:end-1);
    end
    isLog = [false true true true];
%     lbs = [0 0 1 1]; ubs = [0 0 1 1]; ns = 1*ones(1,4);
%     lbs(end) = 5e3; ubs(end) = 5e3; ns(end) = 1; % fix delta_t
    
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
        if fitSpaceOnly
            data.D = data.Ds;
            data.X = sum(data.Xf, 3);
            data.nt = 1;            
        end
        % development/evaluation sets
        [~, evalinds] = reg.trainAndTest(data.X, data.R, devPct);
        [~, foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

        %% run on all cells
        if fitCells
            llstr = 'gauss';
            
            if strcmpi(llstr, 'gauss')                
                lbsc = lbs; ubsc = ubs; nsc = ns; isLogc = isLog;
                scorestr = 'rsq';
            else % remove ssq hyper
                scorestr = 'pseudoRsq';
                ix = [1 3 4]; lbsc = lbs(ix); ubsc = ubs(ix); nsc = ns(ix);
                isLogc = isLog(ix);
            end
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            hyperOpts = struct('lbs', lbsc, 'ubs', ubsc, 'ns', nsc, ...
                'isLog', isLogc);
            opt = struct();
            
            % full ASD and ML
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr);
            eMAP = @(hyper) asd.fitHandle(hyper, data.D, llstr, 'evi', opt);
            Db = asd.sqdist.space(data.Xxy);
            BMAP = @(hyper) asd.fitHandle(hyper, Db, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);

            cell_inds = 1:size(data.Y_all, 2);
%             cell_inds = [5];
            ncells = numel(cell_inds);
            for nn = 1:ncells
                cell_ind = cell_inds(nn);
                if data.neurons{cell_ind}.dPrime < 0.4
                    continue;
                end
                lbl = [data.neurons{cell_ind}.brainArea '-' ...
                    num2str(cell_ind)];

                data.Y = data.Y_all(:,cell_ind); % choose cell for analysis
                if sum(~isnan(data.Y)) == 0
                    continue;
                end
                fits = fitSTRF(data, ML, MAP, BMAP, eMAP, scoreFcn, ...
                    hyperOpts, figdir, lbl, fitMask, ...
                    foldinds, evalinds);

                if ~isempty(fitdir)
                    fits.isLinReg = true;
                    fits.foldinds = foldinds;
                    fits.evalinds = evalinds;
                    fits.dt = datestr(now);
                    tools.updateStruct(dat_fnfcn(lbl), fits);
                end

            end
        end
        %% run on decision
        if fitBehavior
            llstr = 'bern';
            
            if strcmpi(llstr, 'gauss')
                lbsc = lbs; ubsc = ubs; nsc = ns; isLogc = isLog;
                scorestr = 'rsq';
            else % remove ssq hyper
                scorestr = 'pctCorrect'; % pseudoRsq
                ix = [1 3 4];
                if fitSpaceOnly
                    ix = [1 3];
                end
                lbsc = lbs(ix); ubsc = ubs(ix); nsc = ns(ix);
                isLogc = isLog(ix);
            end
            scoreFcn = reg.scoreFcns(scorestr, llstr);
            hyperOpts = struct('lbs', lbsc, 'ubs', ubsc, 'ns', nsc, ...
                'isLog', isLogc);
            
            MAP = @(hyper) asd.fitHandle(hyper, data.D, llstr);
            BMAP = @(hyper) asd.fitHandle(hyper, data.D, llstr, ...
                'bilinear', struct('shape', {{data.ns, data.nt}}));
            ML = @(~) ml.fitHandle(llstr);
            data.Y = data.R;

            lbl = 'decision';
            fits = fitSTRF(data, ML, MAP, BMAP, nan, scoreFcn, hyperOpts, ...
                figdir, lbl, fitMask, foldinds, evalinds);
            if ~isempty(fitdir)
                fits.isLinReg = false;
                fits.foldinds = foldinds;
                fits.evalinds = evalinds;
                fits.dt = datestr(now);
                tools.updateStruct(dat_fnfcn(lbl), fits);
            end
        end
%         close all
    end
    warning on MATLAB:nearlySingularMatrix;
end
