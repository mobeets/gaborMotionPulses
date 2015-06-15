function plotSaccadeKernelOverlay(stim, n, f, showTargs, showHyperflow, lbl)
    if nargin < 6
        lbl = [n.exname ' - ' f.label ' - ' sprintf('%0.2f', f.score)];
    end
    if nargin < 4
        showTargs = true;        
    end
    if nargin < 5
        showHyperflow = true;
    end
    
    
    cmap = cbrewer('div', 'RdBu', 100, 'pchip');
%     cmap = jet(100);
    colormap(cmap);

    xl = nan;
    if isfield(n, 'delayedsaccades') && ~isempty(n.delayedsaccades)
        if strcmp(n.brainArea, 'LIP')
            contourf(n.delayedsaccades.xax, n.delayedsaccades.yax, ...
                n.delayedsaccades.RF, 'LineColor','none');
            caxis([-0.5 1.5]); % less saturated
            xl = [min(n.delayedsaccades.xax) max(n.delayedsaccades.xax)];
            yl = [min(n.delayedsaccades.yax) max(n.delayedsaccades.yax)];
%             imagesc(n.delayedsaccades.xax, n.delayedsaccades.yax, ...
%                 n.delayedsaccades.RF);
        end
    end
    
    set(gca, 'YDir', 'normal'); % imagesc flips y-axis by default
    wf0 = reshape(f.mu(1:end-1), f.shape(1), f.shape(2));
    [wf,~,v] = svds(wf0, 1); % use rank-1 spatial weights
    if sign(wf(1)) ~= sign(v(1))*sign(wf(1))
        wf = -wf;
    end
    
    if n.targPref > 1
        wf = -wf;
        lbl1 = ' (flipped)';
    else
        lbl1 = '';
    end    
    if strcmp(n.brainArea, 'MT')
        sz = 20;
    elseif strcmp(n.brainArea, 'LIP')
        sz = 15;
    else % if contourf shown above
        width = 0.2*norm(median(stim.gaborXY), 2); % in degrees
        sz = sizeToCurUnits(width); % in plot units
    end
    
    t1 = stim.targ1XY;
    t2 = stim.targ2XY;
    if n.targPref > 1
        t3 = t1; t1 = t2; t2 = t3;
    end
    
%     plot.plotKernelSingle(stim.gaborXY, wf(:,1), nan, 3*sz);
%     plot(stim.gaborXY(:,1), stim.gaborXY(:,2), 'ko', 'markersize', sz);
    
    if ~strcmp(n.brainArea, 'LIP') && showHyperflow && ...
            isfield(n, 'hyperflow') && ~isempty(n.hyperflow)
        xl = [min(n.hyperflow.gridx(:)) max(n.hyperflow.gridx(:))];
        yl = [min(n.hyperflow.gridy(:)) max(n.hyperflow.gridy(:))];        
        xs = n.hyperflow.gridx(:);
        ys = n.hyperflow.gridy(:);
        zs1 = n.hyperflow.dx;
        zs2 = n.hyperflow.dy;
        
        proj = @(a,b) dot(a,b)/norm(a);
        dt = diff([median(t1); median(t2)]);
        zs = arrayfun(@(x, y) proj(dt, [x y]), zs1, zs2);
        zs = reshape(zs, numel(unique(xs)), numel(unique(ys)));
        imagesc(unique(xs), unique(ys), -zs);
%         contourf(unique(xs), unique(ys), zs, ...
%             'LineWidth', 1, 'LineColor', 'none');
        caxis([min(0, -max(abs(caxis))) max(0, max(abs(caxis)))]);
        caxis(2*caxis); % less saturated
        
%         keep = arrayfun(@norm, zs1, zs2) >= 2e-2;
%         plot.quiverArrowFix(quiver(xs(keep), ys(keep), ...
%             zs1(keep), zs2(keep)), 80, 120, 'HeadStyle', 'plain', ...
%             'LineWidth', 1.5);

%         quiver(xs(keep), ys(keep), zs1(keep), zs2(keep), ...
%             'k', 'LineWidth', 2);
    end
    
    plot.plotKernelSingle(stim.gaborXY, wf(:,1), nan, 3*sz);
    plot(stim.gaborXY(:,1), stim.gaborXY(:,2), ...
        'ko', 'markersize', sz, 'LineWidth', 2);
    
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
    end
    title([lbl lbl1]);
    
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
