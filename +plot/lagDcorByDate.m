function lagDcorByDate(dt, lags)
% 
% LIP
% 20130515 1 - 1
% 20150313 1 - 3
% 20150326a 2 - 6
% 20150331 1 - 3
% 20150401 2 - 8
% 20150407b 3 - 1
% 
% MT
% 20140304 8 - 1
% 20140305 3 - 2
% 20150304a 6 - 1
% 20150304b 6 - 2
% 20150305b 1 - 3
% 20150306b 2 - 3
% 20150306c 4 - 5
% 20150310 1 - 2
% 20150316c 10 - 3
% 20150324a 6 - 6
% 20150519 1 - 6
% 
% dts = {'20140304', '20140305', '20150304a', '20150304b', '20150305b', '20150306b', '20150306c', '20150310', '20150316c', '20150324a', '20150519'};
% dts = {'20130515',  '20150313', '20150326a', '20150331', '20150401', '20150407b'};
% 
% % dts = {'20150306b', '20150306c', '20150324a', '20150316c'};
% % dts = {'20150326a', '20150401'};
% dts = {'20150304a', '20140307'};
% dts = {'20150304a'};

    isNancy = str2num(dt(4)) >= 5;
    d = io.loadDataByDate2(dt, isNancy);    
    Y0 = d.sps;
    ixType = cellfun(@(n) strcmpi(n.brainArea, 'MT'), d.neurons);
%     ixPref = cellfun(@(n) n.targPref == 1, d.neurons);
%     ixG1 = ixType & ixPref;
%     ixG2 = ixType & ~ixPref;
% 
%     ixG1 = ixType & ~ixPref;
%     ixG2 = ~ixType & ~ixPref;
%     if sum(ixG1) == 0 || sum(ixG2) == 0
%         return;
%     end
    ixG1 = ixType;
    ixG2 = ~ixType;
    A = squeeze(Y0(:,ixG1,:));
    B = squeeze(Y0(:,ixG2,:));    
    [dcs, pvs] = tools.lagDcor(A, B, lags);
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
    title(dt);
    plot(lags, zs0, 'b');
    plot(lags, zs1, 'r');

    figure; hold on;
    plot(lags, nanmean(ws0), 'b');
    plot(lags, nanmean(ws1), 'r');
    title(dt);

%     figure; hold on;
%     ys = dcs{ii};
%     ys(sum(pvs{ii} > 0.05/1)) = nan;
%     ys2 = dcs2{ii};
%     ys2(sum(pvs2{ii} > 0.05/1)) = nan;
%     plot(lags, nanmean(ys));
%     plot(lags, nanmean(ys2), 'r');
%     title(dts{ii});
%     
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
