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
    flipTargPrefNames = {'20150401-LIP_15', ...
        '20150407b-LIP_2', '20150407b-LIP_4', ...
        '20150304a-MT_5', '20150518-MT_5'};
    ignoreMTCells = {'20150407a_25', '20150407a_26', '20150407a_28', ...
        '20150407a_30', '20150407a_31', '20150407a_32', '20140304_14', ...
        '20140307_06', '20140307_08'};
    vals = struct([]);
    for ii = 1:numel(dts)
        dt = dts{ii};
        fs = io.loadFitsByDate(dt, fitdir);
        if isempty(fs)
            continue;
        end
        d = io.loadDataByDate(dt);
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
            if dval.pctCorrect < 0.7
                % skip if monkey not above 70% correct
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
            val.isNancy = str2num(dt(4)) > 4;
            val.dt = dt;
            nm = strsplit(nms{jj}, '_');
            val.type = nm{1};
            if numel(nm) > 1
                val.cellind = str2num(nm{2});
                if any(strcmp([val.dt '_' val.cellind], ignoreMTCells))
                    warning(['Ignoring MT cell ' val.dt '_' val.cellind]);
                    continue;
                end
                val.id = d.neurons{val.cellind}.id;
                val.isLinReg = true;
                val.isCell = true;
                val.llstr = 'gauss';
                val.isLIP = strcmp(val.type, 'LIP');
                val.isMT = strcmp(val.type, 'MT');
%                 continue;
            else
                val.cellind = nan;
                val.id = nan;
                val.isLinReg = false;
                val.isCell = false;
                val.isMT = false;
                val.isLIP = false;
                val.llstr = 'bern';
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
                if isfield(d, 'Y_frz')
                    val.Yfrz = d.Y_frz(:,val.cellind);
                else
                    val.Yfrz = [];
                end
                val.dPrime = neuron.dPrime;
                val.hyper_ssq = val.hyper(2);
                val.targPref = neuron.targPref;
                val.electrodeDepth = neuron.electrodeDepth;
%                 val=getNeuronPrefDir(val, neuron, d.stim);
            else
                val.Y = d.R;
                val.Yfrz = d.R_frz;
                val.targPref = 1;
                val.electrodeDepth = nan;
%                 val=getNeuronPrefDir(val, neuron, d.stim);
            end
            if sum(strcmp(val.name, flipTargPrefNames)) > 0
                val.targPref = ~(val.targPref-1)+1; % flip targ pref
            end
            zerTrials = val.dirprob == 0;
            val.nzertrials = sum(zerTrials);
            val.ntrials = sum(~isnan(val.Y));
            % lowest 20% of full distribution for both monkeys' sessions
            lowTrials = abs(val.dirstrength)/prod(val.shape) < 0.1;%0.015;
            val.nlowmottrials = sum(lowTrials);
            val.fano = nanvar(val.Y)/nanmean(val.Y);
            
            val.C = -d.R+2 == val.targPref;
            val.C0 = d.R; % uncorrected for targPref
            if isfield(d, 'R_frz')
                val.Cfrz = -d.R_frz+2 == val.targPref;
            else
                val.Cfrz = [];
            end
            val.Czer = val.C(zerTrials);
            val.Clow = val.C(lowTrials);
            
            % weights
            if isfield(val, 'w') && numel(val.w) > numel(val.mu)
                val.mu = val.w;
            end
            val.wf = reshape(val.mu(1:end-1), d.ns, d.nt);
            val.b = val.mu(end);
            val = addWeightSubfields(val, val.targPref);
            
            % compare to decision
            if isfield(fs, 'decision')
                fd = fs.decision.(fitstr); fd = fd{end};
                val.wfDec = reshape(fd.mu(1:end-1), d.ns, d.nt);
                wf = val.wf(:);
                if val.targPref == 2
                    wf = -wf;
                end
                val.wfDec_corr = corr(wf, val.wfDec(:));
            end
            
            % responses
            predFcn = reg.getPredictionFcn(val.isLinReg);
            val.Yh = predFcn(d.X, val.mu);
            val.Yh0 = val.Yh - val.b;
            val.Yres = val.Y - val.Yh;
            val.Yzer = val.Y(zerTrials);
            val.Ypos = predFcn(d.X, val.muPos);
            val.Yneg = predFcn(d.X, val.muNeg);
            val.Ylow = val.Y(lowTrials);
            val = tools.autoRegressModelSpikes(val, 4); % creates val.YhAR
            val.YresAR = val.Y - val.YhAR;
            
            % CP
            C = val.C0; % for signed CP (i.e., both above/below 0.5)
            val.cp_Y = tools.AUC(val.Y(C == 1), val.Y(C == 0));
            val.cp_Yres = tools.AUC(val.Yres(C == 1), val.Yres(C == 0));
            val.cp_YresAR = tools.AUC(val.YresAR(C == 1), val.YresAR(C == 0));
            C = val.C;
            val.cp_Yc = tools.AUC(val.Y(C == 1), val.Y(C == 0));
            val.cp_Yresc = tools.AUC(val.Yres(C == 1), val.Yres(C == 0));
            val.cp_YresARc = tools.AUC(val.YresAR(C == 1), val.YresAR(C == 0));
            C = val.Czer;
            val.cp_Yzer = tools.AUC(val.Yzer(C == 1), val.Yzer(C == 0));
            C = val.Clow;
            if any(~isnan(val.Ylow(C == 1))) && any(~isnan(val.Ylow(C == 0)))
                val.cp_Ylow = tools.AUC(val.Ylow(C == 1), val.Ylow(C == 0));
            else
                val.cp_Ylow = nan;
            end
            val = frzTrialWeightedCP(val, d); % check for multiple seeds

            % separability/rank analyses
            if ~isSpaceOnly
                val = addSeparabilityAndRank(val, d);                
            end
            val = centerOfMass(val);
%             if isfield(f, 'hyper')
%                 val = addSelectivityTests(val, f, d, foldinds);
%             end
%             val = compareToPositiveWeights(val, d, foldinds);
            
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
    if isfield(d, 'ix')
        inds = d.ix;
    else
        inds = d.stim.goodtrial & ~d.stim.frozentrials;
    end
    val.ix = inds;
    if isfield(d, 'ixfrz')
        val.ixfrz = d.ixfrz;
    else
        val.ixfrz = [];
    end
    val.dirprob = d.stim.dirprob(inds);
    val.dirstrength = sum(sum(d.stim.pulses(...
        inds,:,:),3),2);
    val.correct = d.stim.correct(inds);
    val.pctCorrect = mean(val.correct);
end

function val = frzTrialWeightedCP(val, d)    
    X = d.stim.pulses(val.ixfrz,:,:);
    Y = val.Yfrz;
    C = val.Cfrz;
    ix = ~isnan(Y); X = X(ix,:,:); Y = Y(ix); C = C(ix);
    val.nfrz = sum(ix);
    
    [~,ix,ic] = unique(mean(X, 3), 'rows');
    cpfrz = zeros(numel(ix),1);
    nf = zeros(numel(ix),1);
    for ii = 1:numel(ix)
        y = Y(ic == ix(ii));
        c = C(ic == ix(ii));
        nf(ii) = sum(ic == ix(ii));
        if sum(c) == 0 || sum(~c) == 0
            continue;
        end
        cpfrz(ii) = tools.AUC(y(c), y(~c));
                
    end
    val.cp_Yfrz = (nf'*cpfrz)/sum(nf);    
    val.nfrz_full = nf';
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

function val = centerOfMass(val)

    % normalized stimulus location
    xs0 = val.Xxy;
    xs = zscore(xs0);

    % thresholded spatial RF
    if isfield(val, 'wfSvd_1')
        ys = val.wfSvd_1(:,1);
    else
        warning('Cannot use wfSvd_1 for centerOfMass. Using w instead.');
        ys = val.w;
    end
    ys = ys*sign(sum(ys)); % make positive
    ys(ys < 0) = 0; % ignore subfields
    ys = ys / sum(ys); % normalize so we can use as weights
    [~,idx] = max(ys);
    
    % mean of location weighted by RF strength
    xc = (xs(:,1)'*ys);
    yc = (xs(:,2)'*ys);
    xc0 = (xs0(:,1)'*ys);
    yc0 = (xs0(:,2)'*ys);
    
    [theta, rho] = cart2pol(xc, yc);
    val.rf_center = [xc yc];
    val.rf_center0 = [xc0 yc0];
    val.rf_center2 = xs(idx,:);
    val.rf_ecc = rho;
    val.rf_theta = theta;
end

function val = addSeparabilityAndRank(val, d)
        
    [U, s, V] = svd(val.wf, 0);
    S = diag(s);
    val.separabilities = S/norm(S,1);
    val.separability_index = val.separabilities(1); % norm-1
    val.inseparability_index = 1 - S(1)/norm(S,2); % norm-2 via Linden
    val.svd_rank = val.shape(2)+1 - sum(cumsum(val.separabilities) >= 0.99);
    predFcn = reg.getPredictionFcn(val.isLinReg);
    
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
    val.YhSvd1 = predFcn(d.X, val.muSvd_1);
    val.YhSvd2 = predFcn(d.X, val.muSvd_2);
    val.YhSvd3 = predFcn(d.X, val.muSvd_3);
    val.Yh0Svd1 = val.YhSvd1 - val.wf(end);
    val.Yh0Svd2 = val.YhSvd2 - val.wf(end);
    val.Yh0Svd3 = val.YhSvd3 - val.wf(end);
end

function val = compareToPositiveWeights(val, d, foldinds)
    if ~strcmp(val.llstr, 'bern')
        return;
    end    
    scoreObj = reg.getScoreObj(false, 'mcc');
    
%     foldinds0 = val.foldinds;
%     X = d.X; Y = d.R; D = d.D;
%     [X, Y, foldinds0] = tools.dropTrialsIfYIsNan(X, Y, foldinds0);
%     obj = reg.getObj_ASD(X, Y, D, scoreObj, struct('foldinds', foldinds0));
%     obj.fitFcnArgFcn = @(obj) {obj.hyper, D, zeros(d.ns*d.nt,1), []};
%     obj.hyperObj = reg.getHyperObj_grid(X, Y, obj, scoreObj);
%     obj = reg.fitAndScore(X, Y, obj, scoreObj);
%     val.wfPos0 = obj.mu;
%     val.wfPos_scores = obj.scores;
%     val.wfPos_hyper = obj.hyper;
    
    foldinds0 = val.foldinds;
    X = d.X; Y = d.R; D = d.D;
    [X, Y, foldinds0] = tools.dropTrialsIfYIsNan(X, Y, foldinds0);
    obj = struct('foldinds', foldinds0, 'hyper', val.hyper);
    obj = reg.getObj_ASD(X, Y, D, scoreObj, obj);
    obj = rmfield(obj, 'hyperObj');
    obj.fitFcnArgFcn = @(obj) {obj.hyper, D, zeros(d.ns*d.nt,1)};
    obj = reg.fitAndScore(X, Y, obj, scoreObj);
    val.wfPos0 = obj.mu;
    val.wfPos0_scores = obj.scores;
    val.wfPos0_hyper = obj.hyper;
    [val.scores; val.wfPos0_scores]
    
%     scsDelta = val.scores - val.wfPos_scores;
%     ms = mean(scsDelta);
%     ses = std(scsDelta)/sqrt(numel(scsDelta));
%     lbs = ms - 2*ses;
%     ubs = ms + 2*ses;
%     val.wfPos_scsDelta = scsDelta;
%     val.wfPos_ms = ms;
%     val.wfPos_lbs = lbs;
%     val.wfPos_ubs = ubs;
%     val.is_better_than_pos = val.wfPos_lbs > 0;
    
    scsDelta = val.scores - val.wfPos0_scores;
    ms = mean(scsDelta);
    ses = std(scsDelta)/sqrt(numel(scsDelta));
    lbs = ms - 2*ses;
    ubs = ms + 2*ses;
    val.wfPos0_scsDelta = scsDelta;
    val.wfPos0_ms = ms;
    val.wfPos0_lbs = lbs;
    val.wfPos0_ubs = ubs;
    val.is_better_than_pos0 = val.wfPos0_lbs > 0;
end

function val = addSelectivityTests(val, f, d, foldinds)    
    
    d.Y = val.Y;
    [ms, lbs, ubs, scs] = tools.rankApprox(f, d, foldinds, val.llstr);
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
    disp(num2str(mean(scs)));

end


function val = getNeuronPrefDir(val, neuron, stim)
% get a neuron's preferred direction
% prefDir = getNeuronPrefDir(neuron)
% Inputs
%   neuron (struct) - output of make_SingleNeuronStruct
% Outputs
%   prefDir - in degrees

val.prefDir   = nan;
val.deltaPref = nan;
val.relativePref = nan;
if isempty(neuron)
    return
end

switch neuron.brainArea
    case 'MT'
        hyperflowPrefDir = nan;
        mtrfmapPrefDir   = nan;
        if ~isempty(neuron.hyperflow) && isfield(neuron.hyperflow, 'prefDir')
            hyperflowPrefDir = neuron.hyperflow.prefDir;
        end
        
        if ~isempty(neuron.mtrfmap)
            mtrfmapPrefDir = neuron.mtrfmap.prefDir;
        end
        
        
        
        val.prefDir = nanmean([hyperflowPrefDir mtrfmapPrefDir]);
        
        
    case 'LIP'
        
        %         if ~isempty(neuron.delayedsaccades)
        %             prefDir = cart2pol(neuron.delayedsaccades.rfXY(1), neuron.delayedsaccades.rfXY(2))*180/pi;
        %         end
        
    otherwise
        warning('getNeuronPrefDir: no support for areas besides MT and LIP')
end



circ_dist = @(x,y) angle(exp(1i*x)./exp(1i*y));
diffCirc  = @(x,y) (circ_dist(x/180*pi, y/180*pi))*180/pi;

val.relativePref=diffCirc(val.prefDir, stim.theta);
val.deltaPref=abs(val.relativePref);
end
