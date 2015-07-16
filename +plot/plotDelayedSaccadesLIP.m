function [xl, yl] = plotDelayedSaccadesLIP(n)
    contourf(n.delayedsaccades.xax, n.delayedsaccades.yax, ...
        n.delayedsaccades.RF, 'LineColor','none');
    caxis([-0.5 1.5]); % less saturated
    xl = [min(n.delayedsaccades.xax) max(n.delayedsaccades.xax)];
    yl = [min(n.delayedsaccades.yax) max(n.delayedsaccades.yax)];
%             imagesc(n.delayedsaccades.xax, n.delayedsaccades.yax, ...
%                 n.delayedsaccades.RF);

end
