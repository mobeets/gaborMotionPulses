n = 21;
C = zeros(n,n);
C(ceil(n/2),:) = linspace(-1,1,n);

fig = figure; hold on;
set(gcf, 'color', 'w');
axis off;
plot.getColors([0 1]);
imagesc(C);

% Create arrow
annotation(fig,'arrow',[0.469367588932806 0.408928571428571],...
    [0.464982778415614 0.292857142857143],'LineWidth',4);

% Create arrow
annotation(fig,'arrow',[0.471428571428571 0.689285714285714],...
    [0.466666666666667 0.633333333333333],'HeadLength',8,'LineWidth',4);

% Create arrow
annotation(fig,'arrow',[0.685633001422475 0.685633001422475],...
    [0.62432915921288 0.486583184257603],'LineStyle',':','LineWidth',2);

% Create arrow
annotation(fig,'arrow',[0.406827880512091 0.408250355618777],...
    [0.298747763864043 0.450805008944544],'LineStyle',':','LineWidth',2);

% Create ellipse
annotation(fig,'ellipse',...
    [0.801853485064011 0.45438282647585 0.0189146514935989 0.0250447227191413],...
    'LineWidth',2,...
    'FaceColor',[0 1 0]);

% Create ellipse
annotation(fig,'ellipse',...
    [0.117728307254624 0.455849731663685 0.0189146514935989 0.0250447227191413],...
    'LineWidth',2,...
    'FaceColor',[0 0.498039215803146 0]);
