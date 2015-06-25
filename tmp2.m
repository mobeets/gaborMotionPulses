vn = tools.makeFitSummaries('data/20150615-nancy/fits', true, 'ASD');
vp = tools.makeFitSummaries('data/20150615-pat/fits', false, 'ASD');
vr = [vp vn];

%%

vt = vr;
fns = fieldnames(vt);
tps = unique({vt.type});
for ii = 1:numel(fns)
    fn = fns{ii};
    if ~isempty(strfind(fn, 'is_'))
        disp(fn);
        for jj = 1:numel(tps)
            vt0 = vt(strcmp({vt.type}, tps{jj}));
            x = numel(vt0([vt0.(fn)]));
            disp([tps{jj} ': ' num2str(x) ' out of ' num2str(numel(vt0))]);
        end
    end
end



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

