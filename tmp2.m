vn2 = tools.makeFitSummaries('data/20150615-nancy/fits', true, 'ASD');
vp2 = tools.makeFitSummaries('data/20150615-pat/fits', false, 'ASD');
vr = [vp2 vn2];

%%

vv = vd;
% vv = vv([vv.isCell]);
vv = vv([vv.is_inseparable0]);
scs = cell2mat({vv.svd_scs}');
% figure; hold on;

xs = scs(:,end);
ys = scs(:,end-2);
zs = [vv.separability_index];
scatter(zs, (ys-xs)./xs, 'b');

%%

scatter(xs, ys);
% xlim([0 1]);
ylim(xlim);
% plot(xlim, ylim, 'k--');
xlabel('ASD'); ylabel('ASD rank-1');


%%

vz = tools.makeFitSummaries('data/20150615-nancy/fits', true, 'ASD', {'20150304a'});

%%

