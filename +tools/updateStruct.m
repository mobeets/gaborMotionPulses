function newObj = updateStruct(filename, obj)
% 
% updates data stored in filename, if any
%   by adding obj to its collection
% 
% EXAMPLE:
% filename contains two previous objs:
%    ML: {2x1 cell}
%    ASD: {2x1 cell}
%    ...
%    test: {2x1 cell}
% 
% while obj looks like:
%    ML: [anything]
%    ASD: [anything]
%    ...
%    test: [anything]
% 
% the resulting obj saved in filename is then:
%    ML: {3x1 cell}
%    ASD: {3x1 cell}
%    ...
%    test: {3x1 cell}
% 
    if exist(filename, 'file') == 2
        oldObj = load(filename);
        newObj = makeUpdatedObj(oldObj, obj);
    else        
        newObj = makeEmptyObj(obj);
    end
    save(filename, '-struct', 'newObj');
end
 
function newObj = makeEmptyObj(obj)
    fieldnms = fieldnames(obj);
    for ii = 1:numel(fieldnms)
        fldnm = fieldnms{ii};
        newObj.(fldnm) = cell(1,1);
        newObj.(fldnm){1} = obj.(fldnm);
    end
end

function newObj = makeUpdatedObj(oldObj, obj)
% 
% updates the fields in oldObj with those in obj,
%   where objs similar to obj are stored in a cell array of oldObj
% 
    fieldnms = fieldnames(oldObj);
    
    % find max index to use
    nobjs0 = nan(numel(fieldnms),1);
    for ii = 1:numel(fieldnms)
        fldnm = fieldnms{ii};
        nobjs0(ii) = numel(oldObj.(fldnm));
    end
    nobjs = max(nobjs0);
    
    % create new struct adding fields in obj to oldObj
    for ii = 1:numel(fieldnms)
        fldnm = fieldnms{ii};
        newObj.(fldnm) = cell(nobjs+1,1); % increase size by one
        curnobjs = numel(oldObj.(fldnm));
        for jj = 1:curnobjs
            newObj.(fldnm){jj} = oldObj.(fldnm){jj};
        end
        item = {};
        if isfield(obj, fldnm)
            item = obj.(fldnm);
        end
        newObj.(fldnm){nobjs+1} = item;
    end
    
    % add fields of obj that aren't already in oldObj
    fieldnms = fieldnames(obj);
    for ii = 1:numel(fieldnms)
        fldnm = fieldnms{ii};
        if ~isfield(newObj, fldnm)
            newObj.(fldnm) = cell(1,1);
            newObj.(fldnm){1} = obj.(fldnm);
        end
    end
end
