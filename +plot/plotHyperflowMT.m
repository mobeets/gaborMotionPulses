function [xl, yl] = plotHyperflowMT(n, t1, t2, contourNoQuiver)
    if n.targPref > 1
        t3 = t2; t2 = t1; t1 = t3;
    end
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
    if contourNoQuiver
        contourf(unique(xs), unique(ys), -zs, ...
        'LineWidth', 1, 'LineColor', 'none');
        caxis([min(0, -max(abs(caxis))) max(0, max(abs(caxis)))]);
        caxis(2*caxis); % less saturated
    else
        imagesc(unique(xs), unique(ys), -zs);
%         keep = arrayfun(@norm, zs1, zs2) >= 2e-2;
        keep = arrayfun(@norm, zs1, zs2) >= 0;
        plot.quiverArrowFix(quiver(xs(keep), ys(keep), ...
            zs1(keep), zs2(keep)), 80, 120, 'HeadStyle', 'plain', ...
            'LineWidth', 1.5);
        quiver(xs(keep), ys(keep), zs1(keep), zs2(keep), ...
            'k', 'LineWidth', 2);
    end
end
