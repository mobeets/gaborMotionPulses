function showTemporalVectors(va, type, nV, addJitter)
    if nargin < 4
        addJitter = false;
    end
    
%     va = va([va.is_selective0]);
    va = va(strcmp({va.type}, type));
    vals = cell2mat(cellfun(@(x) x(:,nV), {va.wfSvd_V}, 'uni', 0))';
    if isempty(vals)
        return;
    end
    scale = sign(vals(:,1));
    if addJitter
        scale = scale + (rand(size(vals,1),1)/1e2);
    end
    vals = vals.*repmat(scale, 1, size(vals,2));
    
    se = std(vals)/sqrt(size(vals,1));

    hold on;
    plot(vals', 'Color', [0.9 0.9 0.9]);
    plot(mean(vals), 'k', 'LineWidth', 3);
    plot.shadedErrorBar([], mean(vals), 2*se, {'LineWidth', 2});    
%     ylim([0 1]);
    title(['Temporal vector #' num2str(nV) ' for ' type ' fits']);
end
