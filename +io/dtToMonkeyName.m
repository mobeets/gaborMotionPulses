function mnkNm = dtToMonkeyName(dt, fullName)
    if nargin < 2
        fullName = true;
    end
    if strcmpi(dt(1:4), '2015')
        mnkNm = 'Nancy';
    else
        mnkNm = 'Pat';
    end
    if ~fullName
        mnkNm = mnkNm(1);
    end
end
