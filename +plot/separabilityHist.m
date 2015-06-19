vv = vf;
tps = {'MT', 'LIP', 'decision'};

figure;
set(gcf, 'DefaultAxesFontSize', 18);
set(gcf, 'color', 'white');
colormap(cbrewer('seq', 'Greens', 2));

for ii = 1:numel(tps)
    vvc = vv(strcmp({vv.type}, tps{ii}));
    subplot(numel(tps),1,ii); hold on;
    hist([vvc.separability_index], linspace(0,1,10));
    xlabel(tps{ii});
    xlim([-0.1 1.1]);
end
subplot(numel(tps),1,1);
title('Separability indices for RF fits');
set(gcf, 'Position', [100 100 400 800]);
plot.saveFigure('separability', 'figures');

%%

figure;
set(gcf, 'DefaultAxesFontSize', 18);
set(gcf, 'color', 'white');

subplot(3,1,1);
plot.showTemporalVectors(vv, 'MT', 1);
ylim([0 1]);
subplot(3,1,2);
plot.showTemporalVectors(vv, 'LIP', 1);
ylim([0 1]);
subplot(3,1,3);
plot.showTemporalVectors(vv, 'decision', 1);
ylim([0 1]);
set(gcf, 'Position', [100 100 400 800]);
plot.saveFigure('temporal-vectors', 'figures');
