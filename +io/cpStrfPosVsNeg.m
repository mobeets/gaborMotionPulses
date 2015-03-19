function vals = cpStrfPosVsNeg(vals)
    inds = [vals.score]./[vals.scoreSdev] > 1;
    inds = [vals.score] > 0.05;
    vals = vals(inds);
    
    [~, inds] = sort({vals.dt});
    vals = vals(inds);
    
    last_dt = '';
    for ii = 1:numel(vals)
        if ~strcmp(last_dt, vals(ii).dt)            
            d = io.loadDataByDate(vals(ii).dt);
            last_dt = vals(ii).dt;
        end
        nm = strsplit(vals(ii).name, '_');
        disp([vals(ii).dt ' ' nm{2}]);
        neuron = d.neurons{str2num(nm{2})};
        
        Xpos = d.X;
        Xpos(Xpos < 0) = 0;
        Xneg = d.X;
        Xneg(Xneg > 0) = 0;
        if neuron.targPref == 2
            Xtmp = Xpos;
            Xpos = Xneg;
            Xneg = Xtmp;
        end
        
        A = getSpikes(Xpos, vals(ii).wf);
        B = getSpikes(Xneg, vals(ii).wf);        
        vals(ii).cp_posneg = tools.AUC(A, B);
        
        C = getSpikes(d.X, vals(ii).wf);
        pref = -(d.R-2) == neuron.targPref;
        C1 = C(pref); % pref choice
        C2 = C(~pref); % anti choice
        vals(ii).cp_pred = tools.AUC(C1, C2);        
        vals(ii).cp = io.getCP(d.stim, neuron, ...
            d.stim.targchosen == nanmax([neuron.targPref, 1]), 0.0, 1.5);
    end
    
    categs = unique({vals.type});
    clrs = lines(numel(categs));
    cp = [vals.cp];
    cp2 = [vals.cp_pred];
%     cp = [vals.cp_pred];
%     cp2 = [vals.cp_posneg];
    figure; hold on;
    xlabel('CP by choice, real spikes');
    ylabel('CP by choice, pred spikes');
%     xlabel('CP by choice');
%     ylabel('CP by (+) vs. (-) stimuli');
    for ii = 1:numel(categs)
        inds = strcmp({vals.type}, categs(ii));
        scatter(cp(inds), cp2(inds), 50, ...
            'DisplayName', categs(ii), 'MarkerFaceColor', clrs(ii,:));
    end
    plot([0 1], [0 1], 'k--');
    legend(categs, 'Location', 'NorthEastOutside');
        
    for ii = 1:numel(categs)
        figure; hold on;
        inds = strcmp({vals.type}, categs(ii));
        hist(cp(inds), linspace(0.0, 1.0, 20));
    end
end

function sps = getSpikes(X, w)
    sps = tools.ifNecessaryAddDC(X, w)*w;
end
