% set xnm, ynm, cnm, fnm

scs = pairs;

% ix = ([scs.cell1_dPrime] > 0.5) & ([scs.cell2_dPrime] > 0.5);
% ix = [scs.rscMatch];
% scs = scs(ix);

samePool = [scs.sameTarg];
% posRfCorr = [scs.sameCorr];
posRfCorr = [scs.sameCorr_spatial];
ix1 = ~samePool & ~posRfCorr;
ix2 = ~samePool & posRfCorr;
ix3 = samePool & ~posRfCorr;
ix4 = samePool & posRfCorr;
ix0 = true(size(ix1));

xnms = {'rfCorr', 'noiseCorrAR', 'rfCorr', 'rfDist'};
ynms = {'noiseCorrAR', 'scoreGainWithCorrs', 'scoreGainWithCorrs', 'noiseCorrAR'};
fnms = {'rscVsRf', 'scoreVsrsc', 'scoreVsRf', 'rfdistVsrsc'};
cnm = 'rfCorr';

% xnms = {'rfDist'};
% ynms = {'noiseCorrAR'};
% fnms = {'rscVsRf'};

fldnms = struct('rfCorr', 'RF correlation', ...
    'rfCorr_spatial', 'RF_{spatial} correlation', ...
    'noiseCorrAR', 'noise correlation', ...
    'rfDist', 'RF center distance', ...
    'scoreGainWithCorrs', '\Delta decoding accuracy');
fldnms = containers.Map(fieldnames(fldnms), struct2cell(fldnms));

mkrsz = 50;
nclrs = 21;
cs1 = cbrewer('seq', 'Reds', nclrs);
cs2 = cbrewer('seq', 'Blues', nclrs);
valToClrInd = @(v, vmn, vmx) floor((nclrs-1)*(v - vmn)/(vmx - vmn))+1;

xlbl = fldnms(xnm);
ylbl = fldnms(ynm);
fnm = fnms{ii};
% fnm = fullfile('data', 'figs', 'summary', [fitnm '_' fnm '.png']);

% figure; hold on; set(gcf, 'color', 'w');
t1 = scs(ix1);
t2 = scs(ix2);
t3 = scs(ix3);
t4 = scs(ix4);
t5 = scs(ix0);
%     c1 = [0.8 0.3 0.3];
%     c2 = [0.3 0.3 0.8];    

if strcmpi(ylbl, '\Delta decoding accuracy')
    scale = 100;
else
    scale = 1;
end

xmx = max(abs([scs.(xnm)]));
ymx = scale*max(abs([scs.(ynm)]));
xl = [-xmx xmx];
yl = [-ymx ymx];

nms = {'different pools', 'same pool', 'both'};
% nrows = 2; ncols = 1;

fig1 = plot.init;
clrs = cs1(valToClrInd([t1.(cnm)], -1, 1),:);
scatter([t1.(xnm)], scale*[t1.(ynm)], mkrsz, clrs, 'filled');
clrs = cs1(valToClrInd([t2.(cnm)], -1, 1),:);
scatter([t2.(xnm)], scale*[t2.(ynm)], mkrsz, clrs, 'filled');
xlabel(xlbl); ylabel(ylbl);
plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
xlim(xl); ylim(yl);
set(gca, 'TickLength', [0 0]);
if ~doSave
    title(nms{1});
else
    title('');
end
plot.setPrintSize(gcf, struct('width', 4, 'height', 3));
if strcmpi(ylbl, '\Delta decoding accuracy')
    ytcks = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(n) [num2str(n) '%'], ytcks, 'uni', 0));
end

fig2 = plot.init;
clrs = cs2(valToClrInd([t3.(cnm)], -1, 1),:);
scatter([t3.(xnm)], scale*[t3.(ynm)], mkrsz, clrs, 'filled');
clrs = cs2(valToClrInd([t4.(cnm)], -1, 1),:);
scatter([t4.(xnm)], scale*[t4.(ynm)], mkrsz, clrs, 'filled');
xlabel(xlbl); ylabel(ylbl);
plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
xlim(xl); ylim(yl);
set(gca, 'TickLength', [0 0]);
if ~doSave
    title(nms{2});
else
    title('');
end
plot.setPrintSize(gcf, struct('width', 4, 'height', 3));
if strcmpi(ylbl, '\Delta decoding accuracy')
    ytcks = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(n) [num2str(n) '%'], ytcks, 'uni', 0));
end

fig3 = plot.init;
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
legend off;
xlabel(xlbl); ylabel(ylbl);
plot(xl, [0 0], '-', 'Color', 0.8*ones(3,1));
plot([0 0], yl, '-', 'Color', 0.8*ones(3,1));
xlim(xl); ylim(yl);
set(gca, 'TickLength', [0 0]);
if ~doSave
    title(nms{3});
else
    title('');
end
plot.setPrintSize(gcf, struct('width', 4, 'height', 3));
if strcmpi(ylbl, '\Delta decoding accuracy')
    ytcks = get(gca, 'YTick');
    set(gca, 'YTickLabel', arrayfun(@(n) [num2str(n) '%'], ytcks, 'uni', 0));
end
