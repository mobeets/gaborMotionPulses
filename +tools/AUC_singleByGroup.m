function [cps, grps] = AUC_singleByGroup(A, B, grps, Agrps, Bgrps)
% 
% calculates CP conditional on group membership
% 
% A, B - [n x nw]
% Agrps, Bgrps - [n x 1]
% grps - [ngrps x 1], or [n x 1] if Agrps/Bgrps not provided
% 
    if nargin < 4
        Agrps = grps;
        Bgrps = grps;
        grps = unique(grps);
    end
    cps = nan(numel(grps), size(A,2));
    for ii = 1:numel(grps)
       grp = grps(ii);
       if sum(Agrps==grp) < 3 || sum(Bgrps==grp) < 3
           continue;
       end
       cps(ii,:) = tools.AUC(A(Agrps==grp,:), B(Bgrps==grp,:));
    end
end
