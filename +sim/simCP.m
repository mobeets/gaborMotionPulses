function obj = simCP(w, xopt, seed)
    if size(w,1) == 1 && size(w,2) > 1
        w = w';
    end
    if nargin > 2
        rng(seed);
    end
    
    % stimulus
    ntrials = 1000;
    ngabors = size(w,1);
    X = makeX(ntrials, ngabors, xopt);
    
    % spikes
    e_sps_chc = 1*randn(ntrials,1); % noise given to both
    % higher gain on e_sps_chc lowers contribution of weights to CP, 
    % and raises contributions of Yres to CP
    e_sps = 2*randn(ntrials,1); % noise given to spikes only
    % higher gain on e_sps shrinks maximum deviations of CP from 0.5 for
    % any signal containing noise, i.e., cp_Yres, cp_Yposres, cp_Ynegres
    Yh = X*w;
    Y = X*w + e_sps + e_sps_chc;    
    wneg = makeWsigned(w, -1);
    Yneg = X*wneg;
    wpos = makeWsigned(w, 1);
    Ypos = X*wpos;
    Yres = e_sps + e_sps_chc;
    Ynegres = Yneg + Yres;
    Yposres = Ypos + Yres;
    
    % choice
    thresh = 0;
    e_chc = 1*randn(ntrials,1); % noise given to choice only
    % higher gain on e_chc shrinks maximum deviations of CP from 0.5 for
    % any signal
    C = makeChoice(sum(X,2), e_sps_chc + e_chc, thresh);
%     C = makeChoice(Y, eps, thresh);

    % CP
    cp.cp_Y = tools.AUC(Y(C), Y(~C));
    cp.cp_Yh = tools.AUC(Yh(C), Yh(~C));
    cp.cp_Yres = tools.AUC(Yres(C), Yres(~C));
    cp.cp_Ypos = tools.AUC(Ypos(C), Ypos(~C));
    cp.cp_Yneg = tools.AUC(Yneg(C), Yneg(~C));
    cp.cp_Ynegres = tools.AUC(Ynegres(C), Ynegres(~C));
    cp.cp_Yposres = tools.AUC(Yposres(C), Yposres(~C));

    % save
    obj.X = X;
    obj.w = w;
    obj.wpos = wpos;
    obj.wneg = wneg;    
    obj.Y = Y;
    obj.C = C;
    obj.Yres = Yres;
    obj.Ypos = Ypos;
    obj.Yneg = Yneg;
    obj.cp = cp;

end

function X = makeX(ntrials, ngabors, xopt, Sigma)
    if nargin < 4 && isnan(xopt)
        Sigma = [1 0.9; 0.9 1];
    else
        Sigma = nan;
    end
    if xopt == 1
        X = repmat(randn(ntrials, 1), 1, ngabors);
    else
        X = randn(ntrials, ngabors);
    end
    if ~any(any(isnan(Sigma))) && ~isempty(Sigma)
        mu = zeros(ngabors,1)';        
        X = repmat(mu, ntrials, 1) + X*chol(Sigma);
    end
end

function wh = makeWsigned(w, sign)
    wh = w;
    if sign < 0
        wh(wh > 0) = 0;
    else
        wh(wh < 0) = 0;
    end
end

function C = makeChoice(Y, eps, thresh)
    C = (Y + eps) > thresh;
end
