function fig = scoreVsHyperparam(hyperind, fitdir)
    dts = io.getDates(fitdir);
    xs = struct();
    for ii = 1:numel(dts)
        fit = io.loadFitsByDate(dts{ii}, fitdir);
        fns = fieldnames(fit);
        for jj = 1:numel(fns)
            fn = fns{jj};
            obj = fit.(fn).ASD{end};
            if numel(obj.hyper) >= hyperind
                val = obj.hyper(hyperind);
                sc = obj.scores(1);
                tp = fn(1:2);
                if isfield(xs, tp)
                    xs.(tp) = [xs.(tp); val sc];
                else
                    xs.(tp) = [val sc];
                end
            end
        end
    end
    fig = figure; hold on;
    xlabel('hyperparam'); ylabel('fit score');
    fns = fieldnames(xs);
    clrs = gray(numel(fns)+1);
    bins = linspace(-1, 1, 20);
    for ii = 1:numel(fns)
        ob = xs.(fns{ii});
        subplot(numel(fns)+1, 1, 1); hold on;
        scatter(ob(:,1), ob(:,2), 'o', ...
            'DisplayName', fns{ii}, ...
            'MarkerFaceColor', clrs(ii,:));
        subplot(numel(fns)+1, 1, ii+1); hold on;
        hist(ob(:,2), bins);
        xlabel('fit score');
        ylabel(fns{ii});
    end
    xscale('log');
    subplot(numel(fns)+1, 1, 1);
    legend('Location', 'northeast');
end
