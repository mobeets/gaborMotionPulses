cs = nan(numel(vg),1);
for ii = 1:numel(vg)
   cs(ii) = nanmean((sign(vg(ii).dirstrength)+1)/2 == vg(ii).C);   
end
%%

bs = linspace(0, 1, 20);

figure; hold on;
for ii = 1:numel(bs)
    boxplot(
    
end
%%

xs = dss;
ys = vars;
ys = [vg.wrat];

ix = ~[vg.isCell];
ix = strcmp({vg.type}, 'LIP');
xs = xs(ix);
ys = ys(ix);

%%
xs = lpp(:,1);
ys = lpp(:,2);

%%
bins = linspace(floor(min(xs)),ceil(max(xs)),10);
[bb,ee]=histc(xs, bins);
% figure;
cents = bins(1:end-1) + diff(bins)/2;
boxplot(ys, ee, 'positions', cents, 'labels', cents);
hold on;
scatter(xs, ys);
% [sxs, sys] = sort(grpstats(ee, ys, @median));
% %        pos(sxs) = 1:6;
% figure;
% boxplot(sxs, sys);
% %        boxplot(xs, ys, 'positions', pos)

%%

XS = [vg.wrat];
YS = [vg.score];

XS = vars;
YS = [vg.cp_Yres];

ix1 = [vg.is_better_than_null_Af];
ix2 = [vg.is_better_than_flat_Af];

tps = {'decision', 'MT', 'LIP'};
figure;
for ii = 1:3
    ix = strcmp({vg.type}, tps{ii});    
    xs = XS(ix);
    ys = YS(ix);    

    subplot(3,1,ii);
    scatter(xs, ys, 'b');
%     xlim([0 1]);
%     ylim([0 1]);
    set(gca, 'FontSize', 14);    
    hold on
    scatter(XS(ix&ix1), YS(ix&ix1), 'b', 'filled')
    scatter(XS(ix&ix2), YS(ix&ix2), 'r', 'filled')
    if ii == 1
%         xlim([0 0.55]);
        legend({'all', 'better than null model', 'better than flat weights'}, 'Location', 'Northeast');        
    end
    xlabel(['proportion of negative weights'' contribution to ' tps{ii}  ' RF']);
    ylabel('fit score');
end

