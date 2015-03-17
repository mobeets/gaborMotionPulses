function [auc, crit] = AUC(A, B)
% [auc, crit] = AUC(A,B)
%
% given spike counts for A and B
% returns the spike count that separates the two groups maximally
% and the area under the roc curve for the two groups
% also: A is target group, B is null group
%
% Jay Hennig (2/1/2011, 3/17/2015)

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
