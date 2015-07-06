function [val, pval] = dcor(x, y, sigTest)
% [val, pval] = dcor(x, y, sigTest)
% 
% calculates distance correlation between x and y
%   returns p-value calculated using either:
%       permutation test: when rows are independent samples
%       shift-permutation test: when rows = time
% 
    if nargin < 3
        sigTest = 'perm';
    end
    df = @(d) d - bsxfun(@plus, mean(d), mean(d,2)) + mean(d(:));    
    A = df(pdist2(x,x));
    B = df(pdist2(y,y));
    cf = @(X, Y) sum(sum(X.*Y))/(size(X,1)^2);
    val = sqrt(cf(A,B)/sqrt(cf(A,A)*cf(B,B)));
    
    if strcmpi(sigTest, 'perm') % permutation test
        N = 1000;
        vals = nan(N,1);
        [~, ix] = sort(rand(N, size(x,1)),2);
        for ii = 1:N
         vals(ii) = tools.dcor(x(ix(ii,:),:), y, false);
        end
        pval = sum(vals > val)/N;
    elseif strcmp(sigTest, 'shift') % shift-permutation test: rows = time
        t0 = ceil(size(x,1)*0.2); % start 20% in
        vals = nan(size(x,1)-t0-1,1);
        for ii = 1:size(x,1)-t0-1
            vals(ii) = tools.dcor(circshift(x, [t0-1+ii 0]), y, false);
        end
        pval = sum(vals > val)/numel(vals);
    end
end
