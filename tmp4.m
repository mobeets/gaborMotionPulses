wf = reshape(v.w, 25, 7);
figure;
for jj = 1:5
    for ii = 1:7
        plot.plotKernelSingle(d.Xxy, wf(:,ii));
        pause(0.02)
    end
%     clf;
%     pause(0.1)
end


%%

vtmp = vr;

vL = vtmp(strcmp({vtmp.type}, 'LIP'));
vrL = cell2mat(cellfun(@(x) var(x(:,1)), {vL.wfSvd_U}, 'uni', 0));
ixL = vrL >= 0.005;
vD = vtmp(strcmp({vtmp.type}, 'decision'));
vrD = cell2mat(cellfun(@(x) var(x(:,1)), {vD.wfSvd_U}, 'uni', 0));
ixD = vrD >= 0.005;

[sum(~ixL) sum(ixL)]
[sum(~ixD) sum(ixD)]

dtD = unique({vD.dt});
vL2 = [];
for ii = 1:numel(dtD)
    vv = vL(ixL & strcmp({vL.dt}, dtD{ii}));
    if numel(vv) > 1
        vL2 = [vL2 vv];
        plot.plotKernelByName(vD(strcmp({vD.dt}, dtD{ii})).name);
        for jj = 1:numel(vv)
            plot.plotKernelByName(vv(jj).name);
        end
    end
end
% '20121107'
% '20130611'
% '20150313'
% '20150331'
% '20150401'
% '20150407b'
% '20150501'

%%

vtmp = vr;
cpn = {'cp_Yres'};%, 'cp_Yzer'};
cpn = {'cp_Yfrz'};
tps = {'MT', 'LIP'};

% vtmp = vtmp([vtmp.score] >= 0.1);

% vars = cell2mat(cellfun(@(x) var(x(:,1)), {vtmp.wfSvd_U}, 'uni', 0));
% vtmp = vtmp(vars >= 0.005);
% [sum(vars < 0.005) numel(vars)]

vgd = vtmp(~[vtmp.isCell]);

% vars = cell2mat(cellfun(@(x) var(x(:,1)), {vgd.wfSvd_U}, 'uni', 0));
% dts = {vgd(vars >= prctile(vars, 33)).dt};
% dts = {vgd([vgd.hyper_delta_space] < 50).dt};
% ix = cell2mat(cellfun(@(dt) sum(strcmp(dts, dt)) > 0, {vtmp.dt}, 'uni', 0));
% vtmp = vtmp(ix);

for jj = 1:numel(cpn)
    figure; clf;    
    
    
    for ii = 1:numel(tps)
        subplot(numel(tps), 1, ii); hold on;
        set(gca, 'FontSize', 14)
        axis equal;
        
        
        vt0 = vtmp(strcmp({vtmp.type}, tps{ii}));
        vt0 = vt0([vt0.cp_Yfrz] > 0 & [vt0.cp_Yfrz] < 1);
        
%         vt0 = vt0([vt0.(cpn{jj})] >= 0.5);
        dtix = cell2mat(cellfun(@(dt) find(strcmp(dt, ...
            unique({vt0.dt}))), {vt0.dt}, 'uni', 0));
%         clrs = cbrewer('qual', 'Set1', max(dtix), 'pchip');
        clrs = gray(max(dtix));
        cps = [vt0.(cpn{jj})];
        
        crs = [vt0.wfDec_corr];
        crs = nan(numel(vt0),1);
        for kk = 1:numel(vt0)            
            vd = vgd(strcmp({vgd.dt}, vt0(kk).dt));
            if numel(vd) == 0
                continue;
            end
            w1 = vt0(kk).w;
            w2 = vd.w;
            if vt0(kk).targPref == 2
                w1 = -w1;
            end
            crs(kk) = corr(w1, w2);
        end
        crs = crs';
        
%         crs = crs(1,2:2:end);
        

        for kk = 1:numel(unique(dtix))
%             scatter(crs(dtix==kk), cps(dtix==kk), 60, clrs(kk,:), 'filled');
            
            xs0 = crs(dtix==kk);
            ys0 = cps(dtix==kk);
            [~,ix0] = sort(xs0);
            plot(xs0(ix0), ys0(ix0), 'color', 'k');%clrs(kk,:));
%             if sum(ix0) > 1
%                 scatter(mean(xs0(ix0)), mean(ys0(ix0)), 50, 'k', 'filled');%clrs(kk,:));
%             end
        end
        ix = ~isnan(crs)&~isnan(cps);

        xs = crs(ix)';
        ys = cps(ix)';
        m = fitlm(xs, ys);
        rsq = m.Rsquared.Ordinary;
        pval = m.Coefficients.pValue(2);
        rho = corr(xs, ys);
        cc = m.coefCI;
        cc = cc(2,:);
        m.plot();
        legend off;
        
        xlabel(['corr(' tps{ii} ' RF, decision RF)']);
        ylabel(cpn{jj});
        
        title([tps{ii} ', Y vs. X: corr = ' num2str(rho) ...
            ', fit r^2 = ' num2str(rsq) ', pval = ' num2str(pval) ...
            ', slope CI = ' num2str(cc)]);
        xlim([-1 1]);
        ylim([0 1]);
    end
end

%%
for kk = 1:numel(vr)
    flps = {'20150401-LIP_15', '20150407b-LIP_2', '20150407b-LIP_4'};
    if sum(strcmp(vr(kk).name, flps)) > 0        
        vr(kk).cp_Yres = 1 - vr(kk).cp_Yres;
        vr(kk).cp_Yfrz = 1 - vr(kk).cp_Yfrz;
        [vr(kk).cp_Yres vr(kk).cp_Yfrz]
    end            
end
%%
vrs = vr([vr.isCell]);
for ii = 1:numel(vrs)
    vd = vr(~[vr.isCell] & strcmp({vr.dt}, vrs(ii).dt));
    w1 = vrs(ii).w;
    if vrs(ii).targPref == 2
        w1 = -w1;
    end
%     [vrs(ii).wfDec_corr corr(w1, vd.w)]
%     vrs(ii).wfDec_corr = corr(w1, vd.w);
    vr(strcmp({vr.name}, vrs(ii).name)).wfDec_corr = corr(w1, vd.w);
end

%%

vtmp = vr;
vars = cell2mat(cellfun(@(x) var(x(:,1)), {vtmp.wfSvd_U}, 'uni', 0));
% vtmp = vtmp(vars >= 0.005);
% vtmp = vtmp([vtmp.hyper_delta_space] < 100);
vtmp = vtmp([vtmp.score] >= 0.1);

cpnm = 'cp_Yfrz';
vtmp2 = vtmp(strcmp({vtmp.type}, 'LIP'));
% vtmp2 = vtmp2([vtmp2.hyper_delta_space] > 100);
vtmp2 = vtmp2([vtmp2.cp_Yfrz] > 0 & [vtmp2.cp_Yfrz] < 1);

figure; hold on;
zs = nan(1,numel(vtmp2));
for ii = 1:numel(vtmp2)
    v = vtmp2(ii);
    xs = v.wfSvd_V(:,1);
    
    vd = vtmp(~[vtmp.isCell] & strcmp({vtmp.dt}, v.dt));
    if numel(vd) == 0
        continue;
    end
    ys = vd.wfSvd_V(:,1);
    
    xs = sign(xs(1))*xs;
    ys = sign(ys(1))*ys;
%     if sum(sign(xs)) ~= sum(sign(ys)) 
    
    if sign(xs(1)) ~= sign(ys(1))
        xs = -xs;
    end
    
    assert(sign(xs(1))==1);
    assert(sign(ys(1))==1);
    zs(ii) = corr(xs, ys);
    scatter(zs(ii), v.(cpnm), 50, 'k', 'filled');

    xlabel('corrcoef wDecision');
%     ylabel(cpn{jj});

%     title([tps{ii} ' ' cpn{jj} ' \rho = ' num2str(rho) ...
%         ' r^2 = ' num2str(rsq) ' pval = ' num2str(pval) ...
%         '   CI = ' num2str(cc)]);
    
end

ws = [vtmp2.(cpnm)];
ix = ~isnan(ws)&~isnan(zs);
xs = zs(ix);
ys = ws(ix);
m = fitlm(zs, ws);
rsq = m.Rsquared.Ordinary;
pval = m.Coefficients.pValue(2);
rho = corr(xs, ys);
cc = m.coefCI;
cc = cc(2,:);

m.Coefficients
m.plot();
legend off;
rsq
cc
pval

%%

