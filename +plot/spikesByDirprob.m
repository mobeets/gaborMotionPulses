function spikesByDirprob(w)
    dps = unique(w.dirprob);
    clrs = jet(numel(dps));
    sz = 50;

    figure; hold on;
    title(w.name);
    for ii = 1:numel(dps)
        inds = (w.dirprob == dps(ii));
        scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:), 'filled');
%         inds = (w.dirprob == dps(ii)) & (w.C == 1);
%         scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:), 'filled');
%         inds = (w.dirprob == dps(ii)) & (w.C == 0);
%         scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:));
    end
    xlabel('Y'); ylabel('Yh'); zlabel('Ypos');
    legend(cellstr(num2str(dps(:))))

    figure; hold on;
    title(w.name);
    for ii = 1:numel(dps)
        inds = (w.dirprob == dps(ii));
        scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:), 'filled');
%         inds = (w.dirprob == dps(ii)) & (w.C == 1);
%         scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:), 'filled');
%         inds = (w.dirprob == dps(ii)) & (w.C == 0);
%         scatter3(w.Y(inds), w.Yh(inds), w.Ypos(inds), sz, clrs(ii,:));
    end
    xlabel('Y'); ylabel('Yh'); zlabel('Yneg');
    legend(cellstr(num2str(dps(:))));
end
