fitdir = '20150716';
fitAllSTRFs(fitdir, false, 'cells ASD');
fitAllSTRFs(fitdir, true, 'cells ASD');
fitAllSTRFs(fitdir, false, 'behavior ASD');
fitAllSTRFs(fitdir, true, 'behavior ASD');
%%
fitdir = '20150715-space-only';
fitAllSTRFs(fitdir, false, 'cells ASD space-only');
fitAllSTRFs(fitdir, true, 'cells ASD space-only');

%%

fitdir = '20150716';
% fitdir = '20150615';
% fitdir = '20150715-space-only';
vn = tools.makeFitSummaries(['data/' fitdir '-nancy/fits'], true, 'ASD');
vp = tools.makeFitSummaries(['data/' fitdir '-pat/fits'], false, 'ASD');
vy = [vp vn];

%%

nms = {vw.name};
vy2 = [];
vw2 = [];
for ii = 1:numel(nms)
    v1 = vw(strcmp({vw.name}, nms{ii}));
    v2 = vx(strcmp({vx.name}, nms{ii}));
    if numel(v1) ~= 1 || numel(v2) ~= 1
        disp(num2str([numel(v1) numel(v2)]));
        continue;
    end
    vw2 = [vw2 v1];
    vy2 = [vy2 v2];
end

%%

% tps = {'is_selective0', 'is_selective1', 'is_selective2', ...
%     'is_inseparable0', 'is_inseparable1', ...
%     'is_better_than_ML', 'is_selective_ML', 'is_selective_subfld0', ...
%     'is_selective_subfld1'};
vtmp = vs;

fns = fieldnames(vtmp);
ix = logical(cell2mat(cellfun(@(x) numel(strfind(x, 'is_'))>0, ...
    fns, 'uni', 0)));
tps = fns(ix);

kds = {'MT', 'LIP', 'decision'};
% kds = {'decision'};
% kds = {'MT', 'LIP'};

disp('===============');
disp('===============');
for ii = 1:numel(tps)    
    if numel(strfind(tps{ii}, '_mean'))>0
        continue;
    end
    disp('------------');
    disp(tps{ii});
    for jj = 1:numel(kds)
        v = vtmp(strcmp({vtmp.type}, kds{jj}));
        x = sum([v.([tps{ii} ''])]);
        disp([kds{jj} ' = ' num2str(x) ' of ' num2str(numel(v))]);
    end
end

%%

v = va(strcmp({va.type}, 'MT'));
A = [v.is_selective];
B = [v.is_inseparable];
% A = [v.nfrozen] > 100;
% B = [v.rank95] > 1;

a=sum(A & B);
b=sum(A & ~B);
c=sum(~A & B);
d=sum(~A & ~B);
disp(['va     :' num2str([a, b, c, d])])


%% 

nm = 'MT';

% if strcmp(nm, 'decision')
%     vtmp = va0;
% else
%     vtmp = va2;
% end
vtmp = va;

% vtmp = vtmp([vtmp.is_inseparable0]);
vtmp = vtmp([vtmp.is_selective0]);
vtmp = vtmp(strcmp({vtmp.type}, nm));
% vtmp = vtmp(cellfun(@(x) size(x,2), {vtmp.sep_ts}) == 1);
ts = cell2mat(cellfun(@(x) x(:,1), {vtmp.wfSvd_V}, 'uni', 0))';
% ws = cell2mat(cellfun(@(x) x(1:19,1), {vtmp.sep_wfs}, 'uni', 0))';
scale = sign(ts(:,1));%./max(abs(ts)')';
% scale = 1./sum(ts,2);
% scale = scale + (rand(size(ts,1),1)/1e2); % jitter
ts = ts.*repmat(scale, 1, size(ts,2));

propstr = 'hyper_delta_time';
vals = [vtmp.(propstr)];
% vals = cell2mat(cellfun(@(x) sum(x(2:end)), {vtmp.separabilities}, 'uni', 0)');
isHigh = vals >= prctile(vals, 50);

figure; colormap gray; hold on; ylim([-1 1])
plot(ts', 'k'); title(nm);
% plot(ts(~isHigh,:)', 'k'); title(nm);
% plot(ts(isHigh,:)', 'r'); title(nm);

%%
figure; hold on;
for ii = 1:size(ts,2)
    subplot(size(ts,2), 1, ii);
    hist(ts(:,ii), linspace(-1, 1, 25));
end
subplot(size(ts,2), 1, 1); title(nm);

%%

for ii = 1:numel(scsP)
    Y1 = scsP(ii).Ys(:,1);
    Y2 = scsP(ii).Ys(:,2);
    X = scsP(ii).stim;
    ix = X > 0;
    [nc, ps] = tools.noiseCorr(Y1, Y2, 30, ix);
    scsP(ii).noiseCorr_pavg = nanmean(ps);
    scsP(ii).noiseCorr_R = nc(1);
    scsP(ii).noiseCorr_L = nc(2);
    scsP(ii).noiseCorr_mx = nanmax(nc);
    scsP(ii).noiseCorr_mn = nanmin(nc);
    scsP(ii).noiseCorr_avg = nanmean(nc);
    
    y1b = nanmean(Y1(ix));
    y2b = nanmean(Y2(ix));
    y1a = nanmean(Y1(~ix));
    y2a = nanmean(Y2(~ix));
    scsP(ii).sigCorr_slope = (y2b-y1b)/(y2a-y1a);
    
    m1 = fitlm(Y1(ix), Y2(ix));
    sl1 = m1.Coefficients.Estimate(2);
    pv1 = m1.Coefficients.pValue(2);
    m2 = fitlm(Y1(~ix), Y2(~ix));
    sl2 = m2.Coefficients.Estimate(2);
    pv2 = m2.Coefficients.pValue(2);
    scsP(ii).noiseCorr_L_slope = sl1;
    scsP(ii).noiseCorr_R_slope = sl2;
    scsP(ii).noiseCorr_L_pval = pv1;
    scsP(ii).noiseCorr_R_pval = pv2;
    scsP(ii).noiseCorr_slope = nanmean([sl1, sl2]);
    
    v = [1 scsP(ii).noiseCorr_slope];
    w = [1 scsP(ii).sigCorr_slope];
    scsP(ii).sigNoiseAngleDeg = acosd(v*w'/(norm(v)*norm(w)));
    scsP(ii).sigNoiseAngleDev = abs(scsP(ii).sigNoiseAngleDeg - 90);
end

