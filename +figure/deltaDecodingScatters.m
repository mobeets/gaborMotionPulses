function fig = deltaDecodingScatters(pairs, xnm, ynm, kind)
% kind [str] refers to decision pool: one of {'different', 'same', 'both'}

scs = pairs;
samePool = [scs.sameTarg];
posRfCorr = [scs.sameCorr_spatial];
ix1 = ~samePool & ~posRfCorr;
ix2 = ~samePool & posRfCorr;
ix3 = samePool & ~posRfCorr;
ix4 = samePool & posRfCorr;
ix0 = true(size(ix1));
t1 = scs(ix1);
t2 = scs(ix2);
t3 = scs(ix3);
t4 = scs(ix4);
t5 = scs(ix0);

cnm = 'rfCorr';

fldnms = struct('rfCorr', 'RF correlation', ...
    'rfCorr_spatial', 'RF_{spatial} correlation', ...
    'noiseCorrAR', 'noise correlation', ...
    'rfDist', 'RF center distance', ...
    'scoreGainWithCorrs', '\Delta decoding accuracy');
fldnms = containers.Map(fieldnames(fldnms), struct2cell(fldnms));
xlbl = fldnms(xnm);
ylbl = fldnms(ynm);

mkrsz = 50;
nclrs = 21;
cs1 = cbrewer('seq', 'Reds', nclrs);
cs2 = cbrewer('seq', 'Blues', nclrs);
valToClrInd = @(v, vmn, vmx) floor((nclrs-1)*(v - vmn)/(vmx - vmn))+1;

if strcmpi(ylbl, '\Delta decoding accuracy')
    scale = 100;
else
    scale = 1;
end

xmx = max(abs([scs.(xnm)]));
ymx = scale*max(abs([scs.(ynm)]));
xl = [-xmx xmx];
yl = [-ymx ymx];

fig = plot.init;

if strcmpi(kind, 'different')
    
    nm = 'different pools';
    clrs = cs1(valToClrInd([t1.(cnm)], -1, 1),:);
    scatter([t1.(xnm)], scale*[t1.(ynm)], mkrsz, clrs, 'filled');
    clrs = cs1(valToClrInd([t2.(cnm)], -1, 1),:);
    scatter([t2.(xnm)], scale*[t2.(ynm)], mkrsz, clrs, 'filled');
    
elseif strcmpi(kind, 'same')
    
    nm = 'same pool';
    clrs = cs2(valToClrInd([t3.(cnm)], -1, 1),:);
    scatter([t3.(xnm)], scale*[t3.(ynm)], mkrsz, clrs, 'filled');
    clrs = cs2(valToClrInd([t4.(cnm)], -1, 1),:);
    scatter([t4.(xnm)], scale*[t4.(ynm)], mkrsz, clrs, 'filled');

elseif strcmpi(kind, 'both')
    
    nm = 'both';
    xs = [t5.(xnm)];
    ys = scale*[t5.(ynm)];
    scatter(xs, ys, mkrsz, 0.5*ones(numel(t5),3), 'filled');
    mdl = fitlm(xs, ys);
    h = mdl.plot;
    h(1).Marker = 'o';
    h(1).Color = 'k';
    h(1).MarkerSize = 4.5;
    h(1).MarkerFaceColor = 0.5*ones(3,1);
    h(1).LineWidth = 1.5;
    h(2).LineWidth = 3;
    h(2).Color = [0.8 0.2 0.2];
    
end

xlabel(xlbl); ylabel(ylbl);
plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
xlim(xl); ylim(yl);
set(gca, 'TickLength', [0 0]);
if ~doSave
    title(nm);
else
    title('');
end
legend off;
if strcmpi(ylbl, '\Delta decoding accuracy')
    ytcks = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(n) [num2str(n) '%'], ytcks, 'uni', 0));
end
plot.setPrintSize(gcf, struct('width', 4, 'height', 3));

end
