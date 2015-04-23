bins = linspace(-1, 1, 30);

figure;

subplot(3, 1, 1); hold on;
xlabel('MT');
xlim([-1, 1]);
hist([vs(strcmp({vs.type}, 'MT')).score], bins);
set(get(gca,'child'),'FaceColor','g','EdgeColor','none');

subplot(3, 1, 2); hold on;
xlabel('LIP');
xlim([-1, 1]);
hist([vs(strcmp({vs.type}, 'LIP')).score], bins);
set(get(gca,'child'),'FaceColor','b','EdgeColor','none');

subplot(3, 1, 3); hold on;
xlabel('decision');
xlim([-1, 1]);
hist([vs(strcmp({vs.type}, 'decision')).score], bins);
set(get(gca,'child'),'FaceColor','r','EdgeColor','none');
