function val = jsDivergence(A, B, flipIfNeg)
% function val = jsDivergence(A, B, flipIfNeg)
% 
% calculates Jensen?Shannon divergence between A and B
%   n.b. 0 <= jsDivergence(A,B) <= 1
% 
% thresholds A and B by ignoring any negative values
% if flipIfNeg, and A or B's sum is negative,
%   the values will be flipped prior to thresholding 
%
    A = A(:);
    B = B(:);
    if any(A<0) || any(B<0)
        if nargin == 3 && flipIfNeg
            A = A*sign(sum(A));
            B = B*sign(sum(B));
        end
        A(A < 0) = 1e-3;
        B(B < 0) = 1e-3;
    end

    A = A/sum(A);
    B = B/sum(B);
    C = (A + B)/2;
    dKL1 = (log2(A) - log2(C))'*A;
    dKL2 = (log2(B) - log2(C))'*B;
    val  = (dKL1 + dKL2 )/2;

end
