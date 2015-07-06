function [dcs, pvs] = lagDcor(d, lags)
% LIP
% 20130515 1 - 1
% 20150313 1 - 3
% 20150326a 2 - 6
% 20150331 1 - 3
% 20150401 2 - 8
% 20150407b 3 - 1
% 
% MT
% 20140304 8 - 1
% 20140305 3 - 2
% 20150304a 6 - 1
% 20150304b 6 - 2
% 20150305b 1 - 3
% 20150306b 2 - 3
% 20150306c 4 - 5
% 20150310 1 - 2
% 20150316c 10 - 3
% 20150324a 6 - 6
% 20150519 1 - 6
% 
    Y0 = d.sps;
%     Y0 = d.sps(sign(sum(d.X,2)) == -1, :,:);
    [nt, ~, ~] = size(Y0);

    ixType = cellfun(@(n) strcmpi(n.brainArea, 'MT'), d.neurons);
    ixPref = cellfun(@(n) n.targPref == 1, d.neurons);
    ixG1 = ixType & ixPref;
    ixG2 = ixType & ~ixPref;
%     ixG1 = ixType;
%     ixG2 = ~ixType;
%     ixG1 = ixType & ~ixPref;
%     ixG2 = ~ixType & ~ixPref;
    if sum(ixG1) == 0 || sum(ixG2) == 0
        return;
    end    
    
    nlags = numel(lags);
    dcs = nan(nt,nlags);
    pvs = nan(nt,nlags);    
    for ii = 1:nt
        if mod(ii, 10) == 0
            disp('.');
        end
        s1 = squeeze(Y0(ii,ixG1,:))';
        s2 = squeeze(Y0(ii,ixG2,:))';        
        ix = ~any(isnan(s1),2) & ~any(isnan(s2),2);
        if sum(ix) == 0
            continue;
        end
        s1 = s1(ix,:);
        s2 = s2(ix,:);
        if all(s1(:) == 0) && all(s2(:) == 0)
            continue;
        end
        for jj = 1:nlags
            lg = lags(jj);
            if lg > 0
                s1t = s1(lg:end,:);
                s2t = s2(1:end-lg+1,:);
            elseif lg == 0
                s1t = s1;
                s2t = s2;
            else
                lg = abs(lg);
                s2t = s2(lg:end,:);
                s1t = s1(1:end-lg+1,:);
            end            
            [dc, pval] = tools.dcor(s1t, s2t, 'shift');
            dcs(ii,jj) = dc;
            pvs(ii,jj) = pval;
        end
    end

end
