fitdir = '20150713';
fitAllSTRFs(fitdir, false, 'cells ASD');
fitAllSTRFs(fitdir, true, 'cells ASD');
fitAllSTRFs(fitdir, false, 'behavior ASD');
fitAllSTRFs(fitdir, true, 'behavior ASD');
%%
fitdir = '20150715-space-only';
fitAllSTRFs(fitdir, false, 'cells ASD space-only');
fitAllSTRFs(fitdir, true, 'cells ASD space-only');

%%

fitdir = '20150713';
% fitdir = '20150615';
% fitdir = '20150715-space-only';
vn = tools.makeFitSummaries(['data/' fitdir '-nancy/fits'], true, 'ASD');
vp = tools.makeFitSummaries(['data/' fitdir '-pat/fits'], false, 'ASD');
vy = [vp vn];

% 20150305a-decision

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
