function cells = makeFitSummaries2(fitdir, dts, fitstr)
    if nargin < 2|| isempty(dts)
        dts = io.getDates(fitdir);
    end
    if nargin < 3
        fitstr = 'ASD';
    end
    flipTargPrefNames = {'20150304a-MT_5', '20150518-MT_5'};
    badCells = {'20150407a_25', '20150407a_26', '20150407a_28', ...
        '20150407a_30', '20150407a_31', '20150407a_32', '20140304_14', ...
        '20140307_6', '20140307_8'};
    
    cells = [];
    for ii = 1:numel(dts)
        dt = dts{ii};
        fs = io.loadFitsByDate(dt, fitdir);
        if isempty(fs)
            continue;
        end
        d = io.loadDataByDate(dt);
        stim = getStimulusInfo(d);
        
        cellnms = fieldnames(fs);
        for jj = 1:numel(cellnms)
            cell = stim;
            assert(~isempty(strfind(cellnms{jj}, 'MT')), 'Found non-MT cell');
            if ~isfield(fs.(cellnms{jj}), fitstr)
                continue;
            end            
            fits = fs.(cellnms{jj}).(fitstr);
            assert(numel(fits) == 1, 'Found multiple fits for cell.');
            curfit = fits{1};
                        
            cell.name = [curfit.dt '-' cellnms{jj}];
            ind = strsplit(cellnms{jj}, '_');
            cell.index = str2num(ind{2});            
            
            cell.mu = curfit.mu;
            cell.w = curfit.w;
            cell.wf = reshape(cell.w, d.ns, d.nt);
            cell.wsep = tools.getSeparableRF(cell.wf);
            cell.rfSpatialVariability = var(cell.wsep.spatial_RF);
            cell.rsq = curfit.score_cvMean;
            cell.rsq_se = curfit.score_cvStd/sqrt(numel(curfit.scores));
            
            neuron = d.neurons{cell.index};
            cell.id = d.neurons{cell.index}.id;
            cell.Y = d.Y_all(:,cell.index);            
            predFcn = reg.getPredictionFcn(curfit.isLinReg);
            cell.Yh = predFcn(d.X, cell.mu);
            
%             cell.dPrime = neuron.dPrime;
%             cell.targPref = neuron.targPref;
            cell.dPrime = tools.dprime(cell.Y, stim.dirstrength > 0);
            cell.targPref = (cell.dPrime < 0) + 1;
            
            cell = tools.RfCenterOfMass(cell); % adds rf_ecc, rf_center, rf_theta
            cell = tools.autoRegressModelSpikes(cell, 4); % adds YhAR
            cell.YresAR = cell.Y - cell.YhAR;
            
            % choice probability
            C = -d.R+2 == cell.targPref;
            cell.cp_YresARc = tools.AUC(cell.YresAR(C == 1), cell.YresAR(C == 0));
            
            nm = [cell.dt '_' num2str(cell.id)];
            if any(ismember(badCells, nm))
                disp(['Skipping bad cell: ' nm]);
                continue;
            end

            cells = [cells cell];
        end
    end
end

function stim = getStimulusInfo(d)
    stim.dt = d.dt;
    stim.X = d.X;
    stim.Xxy = d.Xxy;
    stim.ns = d.ns;
    stim.nt = d.nt;
    stim.ntrials = size(d.X,1);
    
    if isfield(d, 'ix')
        inds = d.ix;
    else
        inds = d.stim.goodtrial & ~d.stim.frozentrials;
    end
    stim.ix = inds;
    stim.dirstrength = sum(sum(d.stim.pulses(...
        inds,:,:),3),2);
    stim.pctCorrect = mean(d.stim.correct(inds));
end
