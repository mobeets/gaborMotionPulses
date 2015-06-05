function [vw, vs, vt] = combineStructArrays(vs, vt)
    fns = setdiff(fieldnames(vs), fieldnames(vt));
    for ii = 1:numel(fns)
       [vt.(fns{ii})] = deal(nan);
    end
    fns = setdiff(fieldnames(vt), fieldnames(vs));
    for ii = 1:numel(fns)
       [vs.(fns{ii})] = deal(nan);
    end
    vw = [vs vt];
end
