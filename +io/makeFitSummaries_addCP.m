function vals = makeFitSummaries_addCP(vals, filterBadScores, filterDecisions)
    if nargin < 2
        filterBadScores = true;
    end
    if nargin < 3
        filterDecisions = true;        
    end
    
    if filterDecisions
        inds = ~strcmp({vals.type}, 'decision');
        vals = vals(inds);
    end
    if filterBadScores
        inds = [vals.score]./[vals.scoreSdev] > 1;
        inds = [vals.score] > 0.05;
        vals = vals(inds);
    end
    
    [~, inds] = sort({vals.dt});
    vals = vals(inds);
    
    last_dt = '';
    for ii = 1:numel(vals)
        if strcmp(vals(ii).type, 'decision')
            continue;
        end
        if ~strcmp(last_dt, vals(ii).dt)            
            d = io.loadDataByDate(vals(ii).dt, vals(ii).isNancy);
            last_dt = vals(ii).dt;
        end
        disp(vals(ii).name);
        neuron = d.neurons{vals(ii).cellind};
        
        % stim strength
        vals(ii).dirprob = d.stim.dirprob(d.stim.goodtrial);
        vals(ii).dirstrength = sum(sum(...
            d.stim.pulses(d.stim.goodtrial,:,:),3),2);
        
        % choice
        C = -d.R+2 == neuron.targPref;
        vals(ii).C = C;
        
        % real spikes
        Y = d.Y_all(:, vals(ii).cellind);
        vals(ii).Y = Y;
        
        % predicted spikes and residual
        Yh = getSpikes(d.X, vals(ii).wf);
        Yres = Y - Yh;
        vals(ii).Yh = Yh;
        vals(ii).Yres = Yres;
        Yh0 = Yh - vals(ii).wf(end);
        vals(ii).rat_Y = max(abs(Yh0))/max(Yh);
        
        % positive/negative parts of RF
        wfPos = vals(ii).wf;
        wfPos(wfPos < 0) = 0;
        wfPos(end) = 0; %vals(ii).wf(end);
        wfNeg = vals(ii).wf;
        wfNeg(wfNeg > 0) = 0;
        wfNeg(end) = 0; %vals(ii).wf(end);
        if neuron.targPref == 2
            tmp = wfPos;
            wfPos = wfNeg;
            wfNeg = tmp;
        end
        vals(ii).wfpos = wfPos;
        vals(ii).wfneg = wfNeg;
        vals(ii).wpos_norm = norm(wfPos,1);
        vals(ii).wneg_norm = norm(wfNeg,1);
        vals(ii).wrat = norm(wfNeg,1)/norm(wfPos,1);
        
        if strcmp(vals(ii).name, '20140304-MT_1')
            x=1;
        end

        % predicted spikes to pos/neg of RF
        Ypos = getSpikes(d.X, wfPos);
        Yneg = getSpikes(d.X, wfNeg);
        Yadj1 = Ypos - Yneg; % monkey should flip wneg response
        Yadj2 = Ypos - Yneg + Yres; % monkey should flip wneg response
        Ypr = Ypos + Yres;
        Ynr = Yneg + Yres;
        vals(ii).Ypos = Ypos;
        vals(ii).Yneg = Yneg;
        
        sinds = ~isnan(Y);
        vals(ii).corr_Yposneg = corr(Ypos, Yneg);
        vals(ii).corr_Yposresneg = corr(Y(sinds)-Ypos(sinds), Yneg(sinds));
        vals(ii).corr_Yposobs = corr(Y(sinds), Ypos(sinds));
        vals(ii).corr_Ynegobs = corr(Y(sinds), Yneg(sinds));
        
        % CP
        vals(ii).cp_Y = tools.AUC(Y(C), Y(~C));
        vals(ii).cp_Yh = tools.AUC(Yh(C), Yh(~C));
        vals(ii).cp_Yres = tools.AUC(Yres(C), Yres(~C));
        vals(ii).cp_Yposres = tools.AUC(Ypr(C), Ypr(~C));
        vals(ii).cp_Ynegres = tools.AUC(Ynr(C), Ynr(~C));
        vals(ii).cp_Ypos = tools.AUC(Ypos(C), Ypos(~C));
        vals(ii).cp_Yneg = tools.AUC(Yneg(C), Yneg(~C));
%         vals(ii).cp_Yadj1 = tools.AUC(Yadj1(C), Yadj1(~C));
%         vals(ii).cp_Yadj2 = tools.AUC(Yadj2(C), Yadj2(~C));

        vals(ii).Y_mean = nanmean(Y);
        vals(ii).Y_var = nanvar(Y);
        vals(ii).Ypos_mean = mean(Ypos);
        vals(ii).Yneg_mean = mean(Yneg);
        vals(ii).Ypos_var = var(Ypos);
        vals(ii).Yneg_var = var(Yneg);
        
        vals(ii).Yposlow = abs(Ypos) < 1;
        vals(ii).Yneglow = abs(Yneg) < 1;
        vals(ii).Yres_var = nanvar(Yres);
        vals(ii).Yposlow_count = sum(vals(ii).Yposlow);
        vals(ii).Yneglow_count = sum(vals(ii).Yneglow);
        vals(ii).Yres_poslow_var = nanvar(Yres(vals(ii).Yposlow));
        vals(ii).Yres_neglow_var = nanvar(Yres(vals(ii).Yneglow));
        
        dps = vals(ii).dirprob;
        udps = unique(dps);
        ys = nan(numel(dps),1);
        for jj = 1:numel(udps)
            y = nanmean(vals(ii).Y(dps == udps(jj)));
            ys(jj) = y;
            vals(ii).(['dps_' num2str((jj))]) = y;
        end
        mdl = fitlm(dps, ys);
        vals(ii).dps_rsq = mdl.Rsquared.Ordinary;
        vals(ii).dps_m = mdl.Coefficients.Estimate(2);
        ix = ~isnan(vals(ii).Y);
        mdl = fitlm(vals(ii).dirprob(ix), vals(ii).Y(ix));
        vals(ii).dp_rsq = mdl.Rsquared.Ordinary;
        vals(ii).dp_m = mdl.Coefficients.Estimate(2);

    end
end

function sps = getSpikes(X, w)
    sps = tools.ifNecessaryAddDC(X, w)*w;
end
