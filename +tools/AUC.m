function cp = AUC(A, B, dim)
% cp = AUC(A, B, DIM)
%
% returns the area under the roc curve for discriminating the two groups
% using a given criterion
% 
% A, B - target and null distributions
% If A, B are matrices:
%       AUC(A, B, DIM) computes AUC along the dimension DIM
%
% Jay Hennig (2/1/2011, 3/17/2015)
% 
    if nargin < 3
        if size(A,1) == 1
            dim = 1;
        else
            dim = 2;
        end
    end
    if dim == 1
        A = A';
        B = B';
    end
    nbins = size(A, 2);
    cp = nan(nbins, 1);
    for ii = 1:nbins
        cp(ii) = AUC_single(A(:,ii), B(:,ii));
    end
end

function [auc, crit] = AUC_single(A, B)

    maxsps = max(max(A), max(B)) + 1; % maximum spike count
    crits = -1:maxsps;
    ncrits = numel(crits);
    
    TPR = zeros(ncrits, 1);
    FPR = zeros(ncrits, 1);
    accur = zeros(ncrits, 1);
    for ii = 1:ncrits
        crit = crits(ii);
        FN = sum(A <= crit); 
        TP = sum(A > crit);
        TN = sum(B <= crit);
        FP = sum(B > crit);
        
        TPR(ii) = TP/(TP+FN); % true positive rate
        FPR(ii) = FP/(FP+TN); % false positive rate
        accur(ii) = (TP+TN)/(TN+TP+FN+FP); % accuracy in decoding
    end

    auc = abs(trapz(FPR,TPR));
    [~, ind] = max(accur);
    crit = crits(ind);
end
