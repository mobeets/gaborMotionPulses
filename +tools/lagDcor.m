function [dcs, pvs] = lagDcor(A, B, lags)
% [dcs, pvs] = lagDcor(A, B, lags)
% 
% computes distance correlations at various lags between 
%   data in A and B
% 
% inputs:
%   A [ntrials x nA x nt]
%   B [ntrials x nB x nt]
%   lags [nlags x 1] - row-shifts of A relative to B
%       where negative lags -> data in A comes prior to data in B
%       and positive lags -> data in B comes prior to data in A
% 
% outputs:
%   dcs - [ntrials x nlags] - distance correlations
%   pvs - [ntrials x nlags] - p-values for distance correlations
%                           from shift permutation test
% 
    nt = size(A,1);
    nlags = numel(lags);
    dcs = nan(nt,nlags);
    pvs = nan(nt,nlags);
    
    for ii = 1:nt
        if mod(ii, 10) == 0
            disp('.');
        end
        s1 = squeeze(A(ii,:,:))';
        s2 = squeeze(B(ii,:,:))';
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
