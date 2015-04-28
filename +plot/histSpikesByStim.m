function histSpikesByStim(X, Y, nbins, xlbl, ylbl)
% example: plot.histSpikesByStim(v.dirprob, v.Y)
% 
    if nargin < 5
        ylbl = 'Y';
    end
    if nargin < 4
        xlbl = 'X';
    end
    if nargin < 3
        nbins = 20;
    end

    figure;
    bins = linspace(min(Y), max(Y), nbins);
    xs = unique(X);

    for ii = 1:numel(xs)
        subplot(numel(xs), 1, ii); hold on;
        xlabel([ylbl ' | ' xlbl ' = ' num2str(xs(ii))]);
        hist(Y(X == xs(ii)), bins);
    end

end
