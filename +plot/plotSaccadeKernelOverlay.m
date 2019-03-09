function plotSaccadeKernelOverlay(stim, n, f, showTargs, showHyperflow, ...
    contourNoQuiver, lbl)
    if nargin < 7
        if ~isstruct(n)
            lbl = [f.label ' = ' sprintf('%0.2f', f.score)];
        else
            lbl = [n.exname ' - ' f.label ' = ' sprintf('%0.2f', f.score)];
        end
    end
    if nargin < 6
        contourNoQuiver = true;
    end
    if nargin < 5
        showHyperflow = true;
    end
    if nargin < 4
        showTargs = true;        
    end
    plot.getColors([0 1]);

    t1 = stim.targ1XY;
    t2 = stim.targ2XY;    
    if ~isstruct(n)
        xl = nan;
    elseif strcmp(n.brainArea, 'LIP') && ...
            isfield(n, 'delayedsaccades') && ~isempty(n.delayedsaccades)
        [xl, yl] = plot.plotDelayedSaccadesLIP(n);
    elseif strcmp(n.brainArea, 'MT') && showHyperflow && ...
            isfield(n, 'hyperflow') && ~isempty(n.hyperflow)
        [xl, yl] = plot.plotHyperflowMT(n, t1, t2, contourNoQuiver);
    else
        xl = nan;
    end
    
    set(gca, 'YDir', 'normal'); % imagesc flips y-axis by default
%     wf0 = reshape(f.mu(1:end-1), f.shape(1), f.shape(2));
%     [wf,~,v] = svds(wf0, 1); % use rank-1 spatial weights
%     if sign(wf(1)) ~= sign(v(1))*sign(wf(1))
%         wf = -wf;
%     end
%     if isstruct(n) && n.targPref > 1
%         wf = -wf;
%         lbl1 = ' (flipped)';
%     else
%         lbl1 = '';
%     end
    
    if ~isstruct(n)
        width = 0.2*norm(median(stim.gaborXY), 2); % in degrees
        sz = sizeToCurUnits(width); % in plot units
        sz = 20;
    elseif strcmp(n.brainArea, 'MT')
        sz = 20;
    elseif strcmp(n.brainArea, 'LIP')
        sz = 15;
    end
    
    wf = f.wsep.spatial_RF;
    if n.targPref > 1
        wf = -wf;
    end
    plot.plotKernelSingle(stim.gaborXY, wf, nan, 3*sz);
    plot(stim.gaborXY(:,1), stim.gaborXY(:,2), ...
        'ko', 'markersize', sz, 'LineWidth', 2);
    
    if isstruct(n) && n.targPref > 1
        t3 = t1; t1 = t2; t2 = t3;
    end
    if showTargs        
        if var(t1,1) == 0
            t1 = t1 + 0.1*randn(size(t1));
            t2 = t2 + 0.1*randn(size(t2));
        end
        if ~isnan(xl)
            xl = [min([xl, t1(:,1)', t2(:,1)']) max([xl, t1(:,1)', t2(:,1)'])];
            yl = [min([yl, t1(:,2)', t2(:,2)']) max([yl, t1(:,2)', t2(:,2)'])];
        end
        plot(t1(:,1), t1(:,2), 'o', 'color', [0.2 0.8 0.2]);
        plot(t2(:,1), t2(:,2), 'o', 'color', [0.2 0.5 0.2]);
        scatter(0, 0, 50, [0.8 0.2 0.2], 'filled');
    end
    if ~isnan(xl)
        xlim(xl); ylim(yl);
    end
end

function sz = sizeToCurUnits(sz0)
    curUnits = get(gca, 'units');
    set(gca, 'units', 'points');
    axpos = get(gca, 'Position');
    sz = sz0/diff(xlim)*axpos(3);
    set(gca, 'units', curUnits);
end
