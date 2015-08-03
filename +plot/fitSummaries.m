function S = fitSummaries(vals, x, y, z, plotUnity)
    if nargin < 5
        plotUnity = false;
    end
    if nargin < 4 || (numel(z) == 1 && isnan(z))
        z = x;
    end
    categ = {vals.type};
    categs = unique(categ);
    monkey = [vals.isNancy];
    monkeys = unique(monkey);
    
    clrs = lines(numel(categs)+2);
    hold on;
    
    xs = [vals.(x)];
    ys = [vals.(y)];
    zs = [vals.(z)];
    S = struct();
    S.(x) = xs;
    S.(y) = ys;
    S.(z) = zs;
    S.monkey = monkey;
    S.categ  = categ;
    
    for jj = 1:numel(monkeys)
        mnk_ind = monkey == monkeys(jj);
        if jj == 1
            clrtyp = 'o';
        else
            clrtyp = 's';
        end
        for ii = 1:numel(categs)
            if numel(vals) > 1
                ind = strcmp(categ, categs(ii));
            else
                ind = true(numel(xs),1);
            end
            ind = ind & mnk_ind;
            scatter3(xs(ind), ys(ind), zs(ind), 10, clrs(ii,:), ...
                'Marker', clrtyp, 'DisplayName', categs{ii});
        end
    end
    lbl = @(x) strrep(x, '_', '-');
    xlabel(lbl(x)); ylabel(lbl(y)); zlabel(lbl(z));
    legend(categs, 'Location', 'Best');
    if plotUnity
        plot([0 1], [0 1], 'k--');
    end
end
