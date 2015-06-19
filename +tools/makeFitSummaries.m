function vals = makeFitSummaries(fitdir, isNancy, fitstr, dts)
    if nargin < 1
        fitdir = 'fits';
    end
    if nargin < 2
        isNancy = false;
    end
    if nargin < 3
        fitstr = 'ASD';
    end
    if nargin < 4 || isempty(dts)
        dts = io.getDates(fitdir);
    end
    isSpaceOnly = false;
    knownProblemFields = {'score_noCV', 'muCorrFolds', 'mus'};
    
    vals = struct([]);
    for ii = 1:numel(dts)
        dt = dts{ii};
        fs = io.loadFitsByDate(dt, fitdir);
        if isempty(fs)
            continue;
        end
        d = io.loadDataByDate(dt, isNancy);
        if isSpaceOnly
            d.X = sum(d.Xf, 3);
            d.D = d.Ds;
            d.nt = 1;
        end
        dval = addStimulusInfo(d);
        
        [~, foldinds] = tools.trainAndTestKFolds(d.X, d.R, 10);
        nms = fieldnames(fs);
        for jj = 1:numel(nms)
            if ~isfield(fs.(nms{jj}), fitstr)
                continue;
            end
            
            % load fit
            fs0 = fs.(nms{jj}).(fitstr);
            fs0 = fs0(~cellfun(@isempty, fs0));
            f = fs0{end};
            clear val;
            
            % add fit and stim info
            val = struct();
            val = fillStructWithOtherStruct(val, f);
            val = fillStructWithOtherStruct(val, dval);
            
            % meta data
            val.name = [dt '-' nms{jj}];
            disp(val.name);
            val.isNancy = isNancy;
            val.dt = dt;
            nm = strsplit(nms{jj}, '_');
            val.type = nm{1};
            if numel(nm) > 1
                val.cellind = str2num(nm{2});
                val.isLinReg = true;
                val.isCell = true;
                val.llstr = 'gauss';
%                 continue;
            else
                val.cellind = nan;
                val.isLinReg = false;
                val.isCell = false;
                val.llstr = 'bern';
%                 continue;
            end
            
            % hypers
            if ~isfield(val, 'hyper') || all(isnan(val.hyper))
                val.hyper = nan(3,1);
            end
            val.dPrime = nan;
            val.hyper_ssq = nan;
            val.hyper_ro = val.hyper(1);
            val.hyper_delta_space = val.hyper(end-1);
            val.hyper_delta_time = val.hyper(end);
            
            % scores
            val.score_mean = mean(val.scores);

            % trials
            if val.isCell
                neuron = d.neurons{val.cellind};
                val.Y = d.Y_all(:,val.cellind);
                val.dPrime = neuron.dPrime;
                val.hyper_ssq = val.hyper(2);
                targPref = neuron.targPref;
            else
                val.Y = d.R;
                targPref = 1;                
            end
            val.C = logical(-d.R+2 == targPref);
            val.ntrials = sum(~isnan(val.Y));
            val.fano = nanvar(val.Y)/nanmean(val.Y);
            
            % weights
            if isfield(val, 'w') && numel(val.w) > numel(val.mu)
                val.mu = val.w;
            end
            val.wf = reshape(val.mu(1:end-1), d.ns, d.nt);
            val.b = val.mu(end);
            val = addWeightSubfields(val, targPref);
            
            % responses
            val.Yh = getResponse(d.X, val.mu, val.isLinReg);
            val.Yh0 = val.Yh - val.b;
            val.Yres = val.Y - val.Yh;
            val.Yzer = val.Y(val.dirprob == 0);
            val.Ypos = getResponse(d.X, val.muPos, val.isLinReg);
            val.Yneg = getResponse(d.X, val.muNeg, val.isLinReg);
            
            % CP
            val.cp_Y = tools.AUC(val.Y(val.C), val.Y(~val.C));
            val.cp_Yzer = tools.AUC(val.Yzer(val.C(val.dirprob == 0)), ...
                val.Yzer(~val.C(val.dirprob == 0)));
            val.cp_Yres = tools.AUC(val.Yres(val.C), val.Yres(~val.C));
            
            % separability/rank analyses
            if ~isSpaceOnly
                val = addSeparabilityAndRank(val, d);
            end
            if isfield(f, 'hyper')
                val = addSelectivityTests(val, f, d, foldinds);
            end
            
            % add nans to any fields not present, ignore the rest
            for kk = 1:numel(knownProblemFields)
                if ~isfield(val, knownProblemFields{kk})
                    val.(knownProblemFields{kk}) = nan;
                end
            end
            if ~isempty(vals)
                [vals, val] = dropNonMatchingFields(vals, val);
            end
            vals = [vals val];
        end
    end
end

function x = fillStructWithOtherStruct(x, y)
    fns = fieldnames(y);
    for ii = 1:numel(fns)
        if ~isfield(x, fns{ii})
            x.(fns{ii}) = y.(fns{ii});
        end
    end
end

function [v1, v2] = dropNonMatchingFields(v1, v2)
    fn1 = fieldnames(v1);
    fn2 = fieldnames(v2);
    fn1_rm = setdiff(fn1, fn2);
    fn2_rm = setdiff(fn2, fn1);
    if ~isempty([fn1_rm fn2_rm'])
        disp(['Ignoring fields = ' strjoin([fn1_rm fn2_rm'], ', ')]);
    end
    for ii = 1:numel(fn1_rm)
        v1 = rmfield(v1, fn1_rm{ii});
    end
    for ii = 1:numel(fn2_rm)
        v2 = rmfield(v2, fn2_rm{ii});
    end
end

function val = addStimulusInfo(d)
    val.Xxy = d.Xxy;
    inds = d.stim.goodtrial & ~d.stim.frozentrials;
    val.dirprob = d.stim.dirprob(inds);
    val.dirstrength = sum(sum(d.stim.pulses(...
        inds,:,:),3),2);
end

function Yh = getResponse(X, w, isLinReg)
    Yh = tools.ifNecessaryAddDC(X, w)*w;
    if ~isLinReg
        Yh = tools.logistic(Yh);
    end
end

function val = addWeightSubfields(val, targPref)
    muPos = val.mu; muPos(muPos < 0) = 0; muPos(end) = 0;
    muNeg = val.mu; muNeg(muNeg > 0) = 0; muNeg(end) = 0;
    if targPref == 2
        tmp = muPos; muPos = muNeg; muNeg = tmp;
    end
    val.muPos = muPos;
    val.muNeg = muNeg;
    val.muPos_norm = norm(muPos,1);
    val.muNeg_norm = norm(muNeg,1);
    val.wrat = val.muNeg_norm / (val.muNeg_norm + val.muPos_norm);
end

function val = addSeparabilityAndRank(val, d)
    
    [U, s, V] = svd(val.wf, 0);
    S = diag(s);
    val.separabilities = S/norm(S,1);
    val.separability_index = val.separabilities(1); % norm-1
    val.inseparability_index = 1 - S(1)/norm(S,2); % norm-2 via Linden
    val.svd_rank = val.shape(2)+1 - sum(cumsum(val.separabilities) >= 0.99);
    
    wf_svdf = @(jj) U(:,jj)*S(jj)*V(:,jj)';
    mu_svdf = @(jj) [reshape(wf_svdf(jj), [], 1); val.mu(end)];
    val.wfSvd_U = U;
    val.wfSvd_V = V;
    val.wfSvd_S = S;
    val.wfSvd_1 = wf_svdf(1);
    val.wfSvd_2 = wf_svdf(2);
    val.wfSvd_3 = wf_svdf(3);
    val.muSvd_1 = mu_svdf(1);
    val.muSvd_2 = mu_svdf(2);
    val.muSvd_3 = mu_svdf(3);
    val.YhSvd1 = getResponse(d.X, val.muSvd_1, val.isLinReg);
    val.YhSvd2 = getResponse(d.X, val.muSvd_2, val.isLinReg);
    val.YhSvd3 = getResponse(d.X, val.muSvd_3, val.isLinReg);
    val.Yh0Svd1 = val.YhSvd1 - val.wf(end);
    val.Yh0Svd2 = val.YhSvd2 - val.wf(end);
    val.Yh0Svd3 = val.YhSvd3 - val.wf(end);
end

function val = addSelectivityTests(val, f, d, foldinds)

    d.Y = val.Y;
    [ms, lbs, ubs, scs] = tools.rankApprox2(f, d, foldinds, val.llstr);
    val.svd_ms = ms;
    val.svd_lbs = lbs;
    val.svd_ubs = ubs;
    val.svd_scs = mean(scs);

%     scsDelta = [scsFlatDelta scsNullDelta ...
%         scsML-scsFlat scsML-scsNull ...
%         scs(:,3)-scs(:,1) scs(:,2)-scs(:,1) scs(:,3)-scsML];

    names = {'is_better_than_flat_A1', ...
        'is_better_than_flat_A3', ...
        'is_better_than_flat_Af', ...
        'is_better_than_null_A1', ...
        'is_better_than_null_A3', ...
        'is_better_than_null_Af', ...
        'is_better_than_flat_ML', ...
        'is_better_than_null_ML', ...
        'is_inseparable_Af', ...
        'is_inseparable_A3', ...
        'is_better_than_ML_Af'};
    
    tests = ubs < -1e-3;
    assert(numel(names) == numel(tests));
    for ii = 1:numel(names)
        val.(names{ii}) = tests(ii);
    end
    
    tests0 = ms < -1e-3;
    assert(numel(names) == numel(tests0));
    for ii = 1:numel(names)
        val.([names{ii} '_mean']) = tests0(ii);
    end
    
%     scsDelta = [scs(:,2)-scs(:,3) scs(:,2)-scs(:,4) scs(:,1)-scs(:,4) ...
%             scsDelta scs(:,4)-scsML' scsML'-nullSc];

%     val.is_inseparable0 = ubs(1) < z; % ASD rank-3 better than ASD rank-1
%     val.is_inseparable1 = ubs(2) < z; % ASD        better than ASD rank-1
%     val.is_selective_subfld0 = ubs(3) < z; % thresh(ASD) better than ASD
%     val.is_selective_subfld1 = ubs(4) < z; % thresh(ASD) better than null model
%     val.is_selective0 = ubs(5) < z; % ASD rank-1    better than null model
%     val.is_selective1 = ubs(6) < z; % ASD rank-3    better than null model
%     val.is_selective2 = ubs(7) < z; % ASD           better than null model
%     val.is_better_than_ML = ubs(8) < z; % ASD       better than ML
%     val.is_selective_ML = ubs(9) < z; % ML          better than null model
%     val.is_selective_ML1 = ubs(10) < z; % ML        better than null model
%     
%     assert(numel(ms)==10);
%     z = -1e-3;
%     val.is_inseparable0 = lbs(1) > -z; % ASD rank-3 better than ASD rank-1
%     val.is_inseparable1 = lbs(2) > -z; % ASD        better than ASD rank-1
%     val.is_selective_subfld0 = ubs(3) < z; % thresh(ASD) better than ASD
%     val.is_selective_subfld1 = ubs(4) < z; % thresh(ASD) better than null model
%     val.is_selective0 = ubs(5) < z; % ASD rank-1    better than null model
%     val.is_selective1 = ubs(6) < z; % ASD rank-3    better than null model
%     val.is_selective2 = ubs(7) < z; % ASD           better than null model
%     val.is_better_than_ML = ubs(8) < z; % ASD       better than ML
%     val.is_selective_ML = ubs(9) < z; % ML          better than null model
%     val.is_selective_ML1 = ubs(10) < z; % ML        better than null model
%         
%     val.is_inseparable0a = ms(1) > -z;
%     val.is_inseparable1a = ms(2) > -z;
%     val.is_selective_subfld0a = ms(3) < z;
%     val.is_selective_subfld1a = ms(4) < z;
%     val.is_selective0a = ms(5) < z;
%     val.is_selective1a = ms(6) < z;
%     val.is_selective2a = ms(7) < z;
%     val.is_better_than_MLa = ms(8) < z;
%     val.is_selective_MLa = ms(9) < z;
%     
    disp(num2str(mean(scs)));

end
