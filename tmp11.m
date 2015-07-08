dts = unique({vs.dt});
for ii = 1:numel(dts)
    v0 = vs(strcmp({vs.dt}, dts{ii}) & strcmp({vs.type}, 'LIP'));
    v1 = vs(strcmp({vs.dt}, dts{ii}) & strcmp({vs.type}, 'MT'));
    
    if numel(v0) > 0 && numel(v1) > 0
        disp([dts{ii} ' ' num2str(numel(v0)) ' - ' num2str(numel(v1))]);
    end
    
%     a = sum([v0.targPref] == 1);
%     if a > 0 && numel(v0) > a
%         disp([dts{ii} ' ' num2str(a) ' - ' num2str(numel(v0)-a)]);    
%     end
    
end

%%

dts = {'20140304', '20140305', '20150304a', '20150304b', '20150305b', '20150306b', '20150306c', '20150310', '20150316c', '20150324a', '20150519'};
dts = {'20130515',  '20150313', '20150326a', '20150331', '20150401', '20150407b'};

% dts = {'20150306b', '20150306c', '20150324a', '20150316c'};
% dts = {'20150326a', '20150401'};
dts = {'20150304a', '20140307'};
dts = {'20150304a'};

lags = 0:20;
dcs = cell(numel(dts),1);
pvs = cell(numel(dts),1);
for ii = 1:numel(dts)
    d = io.loadDataByDate2(dts{ii}, str2num(dts{ii}(4)) >= 5);
    [dcs{ii}, pvs{ii}] = tools.lagDcor(d, lags);
    ix = sign(sum(d.X,2)) == -1;
    ix = d.R > -1;
    
    xs = dcs{ii};
    xs0 = xs(ix,:);
    xs1 = xs(~ix,:);

    ys = pvs{ii};
    ys0 = ys(ix,:);
    ys1 = ys(~ix,:);
    
    ix0 = ys0 < 0.05/numel(lags);
    ix1 = ys1 < 0.05/numel(lags);
    
    zs0 = sum(ix0)/size(ys0,1);
    zs1 = sum(ix1)/size(ys1,1);
    ws0 = xs0;
    ws1 = xs1;
    ws0(~ix0) = nan;
    ws1(~ix1) = nan;
    
    figure; hold on;
    title(dts{ii});
    plot(lags, zs0, 'b');
    plot(lags, zs1, 'r');
    
    figure; hold on;
    plot(lags, nanmean(ws0), 'b');
    plot(lags, nanmean(ws1), 'r');
    title(dts{ii});
    
%     figure; hold on;
%     ys = dcs{ii};
%     ys(sum(pvs{ii} > 0.05/1)) = nan;
%     ys2 = dcs2{ii};
%     ys2(sum(pvs2{ii} > 0.05/1)) = nan;
%     plot(lags, nanmean(ys));
%     plot(lags, nanmean(ys2), 'r');
%     title(dts{ii});
    
% 
%     figure; hold on;
%     ys = sum(pvs2{ii} < 0.05/1)/size(pvs2{ii},1);
%     plot(lags, ys);
%     title(dts{ii});
%     
%     figure; hold on;
%     ys = dcs2{ii};
%     ys(sum(pvs2{ii} > 0.05/1)) = 0;
%     plot(lags, nanmean(ys));
%     ys2 = dcs{ii};
%     ys2(sum(pvs{ii} > 0.05/1)) = 0;
%     plot(lags, nanmean(ys) - nanmean(ys2));
%     title(dts{ii});
end
