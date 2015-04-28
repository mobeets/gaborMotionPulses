diststr = 'separability';
xmin = 0;
xmax = 1;
bins = linspace(xmin, xmax, 20);

figure;

subplot(3, 1, 1); hold on;
xlabel('MT');
xlim([xmin, xmax+0.1]);
hist([vs(strcmp({vs.type}, 'MT')).(diststr)], bins);
set(get(gca,'child'),'FaceColor','g','EdgeColor','none');

subplot(3, 1, 2); hold on;
xlabel('LIP');
xlim([xmin, xmax+0.1]);
hist([vs(strcmp({vs.type}, 'LIP')).(diststr)], bins);
set(get(gca,'child'),'FaceColor','b','EdgeColor','none');

subplot(3, 1, 3); hold on;
xlabel('decision');
xlim([xmin, xmax+0.1]);
hist([vs(strcmp({vs.type}, 'decision')).(diststr)], bins);
set(get(gca,'child'),'FaceColor','r','EdgeColor','none');
