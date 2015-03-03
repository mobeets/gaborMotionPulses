function f = colorScheme(clrNeg, clrMid, clrPos)
    if nargin < 1
        clrPos = [0.3, 0.3, 0.9];
        clrNeg = [0.9, 0.3, 0.3];
        clrMid = [0.95, 0.95, 0.95];
    end
    f = @(x) getColor(x, clrNeg, clrMid, clrPos);
end

function v = getColor(x, clrNeg, clrMid, clrPos)
    if x >= 0
        v = getColor2(x, clrMid, clrPos);
    else
        v = getColor2(-x, clrMid, clrNeg);
    end 
end

function v = getColor2(x, clrStart, clrEnd)
    v = clrStart + x*(clrEnd - clrStart);
end
