function c = cellColsTocellStr(a, b, delim)
    ab = [asStr(a); asStr(b)]';
    f = @(x) strjoin(x, delim);
    c = cellfun(f, mat2cell(ab, ones(size(ab,1),1), 2)', 'uni', 0);
end

function y = asStr(x)
% converts x to cell array of str, if necessary
    y = x;
    if islogical(x) || islogical(x)
        y = cellstr(num2str(x(:)))';
    end
end
