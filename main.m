%% load
data = io.loadData('+io/XY.mat');
cell_names = [14 15 16 19 20 21 23 24 25 27];
% i.e. column i of data.Y ~> cell # cell_names(i)

%% params

nfolds = 5;
fold_for_plots = 1;
[~,~,~,~,foldinds] = reg.trainAndTestKFolds(data.X, data.R, nfolds);

figdir = 'figs';
datdir = 'fits';
dat_fnfcn = @(tag) fullfile(datdir, [tag '.mat']);
fig_fnfcn = @(tag, ext) fullfile(figdir, [tag '.' ext]);
fig_svfcn = @(fig, tag, ext) hgexport(fig, fig_fnfcn(tag, ext), hgexport('factorystyle'), 'Format', ext);

%% run on all cells

isLog = true;
isLinReg = true;
llstr = 'gauss';
% for gridding:
lbs = [-3, -5 -5]; ubs = [3, 10 10]; ns = 5*ones(1,3);
hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
M = asd.linearASDStruct(data.D, llstr);
mlFcn = @(~) ml.fitopts('gauss'); % no poisson for ML yet

% cell_inds = 2;
cell_inds = 1:size(data.Y_all, 2);

ncells = numel(cell_inds);
for nn = 1:ncells
    cell_ind = cell_inds(nn);
    lbl = ['cell_' num2str(cell_ind)];
    data.Y = data.Y_all(:,cell_ind); % choose cell for analysis

    [fits.ASD, wf, sc] = reg.cvMaxScoreGrid(data, hypergrid, M.mapFcn, {}, ...
        M.rsqFcn, {}, foldinds, ['ASD-' llstr], fold_for_plots);
    fits.ASD.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
        data.nt, fold_for_plots, fits.ASD.lbl, sc);
    fig_svfcn(fits.ASD.fig, [lbl '-ASD_' llstr], 'png');

%     [fits.ASD_gs, wf, sc] = reg.cvMaxScoreGridSearch(data, lbs, ubs, ns, M.mapFcn, ...
%         {}, M.rsqFcn, {}, foldinds, fold_for_plots, 'ASD-gs', isLog);
%     fits.ASD_gs.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
%         data.nt, fold_for_plots, fits.ASD_gs.lbl, sc);
%     fig_svfcn(fits.ASD_gs.fig, [lbl '-ASD-gs_' llstr], 'png');

    [fits.ML, wf, sc] = reg.cvMaxScoreGrid(data, [nan nan nan], mlFcn, {}, ...
        M.rsqFcn, {}, foldinds, 'ML', 1);
    fits.ML.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
        data.nt, fold_for_plots, fits.ML.lbl, sc);
    fig_svfcn(fits.ML.fig, [lbl '-ML'], 'png');

    fits.isLinReg = isLinReg;
    io.updateStruct(dat_fnfcn(lbl), fits);
        
end

%% run on decision

isLog = true;
isLinReg = false;
lbs = [-3, -2, -5 -5]; ubs = [3, 10, 10 10]; ns = 5*ones(1,4);
hypergrid = exp(tools.gridCartesianProduct(lbs, ubs, ns));
M = asd.logisticASDStruct(data.D);
mlFcn = @(~) ml.fitopts('bern');
data.Y = data.R;

[fits.ASD, wf, sc] = reg.cvMaxScoreGrid(data, hypergrid, M.mapFcn, {}, ...
    M.rsqFcn, {}, foldinds, 'ASD', fold_for_plots);
fits.ASD.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
    data.nt, fold_for_plots, fits.ASD.lbl, sc);
fig_svfcn(fits.ASD.fig, 'decision-ASD', 'png');

% [fits.ASD_gs, wf, sc] = reg.cvMaxScoreGridSearch(data, lbs, ubs, ns, M.mapFcn, ...
%     {}, M.rsqFcn, {}, foldinds, fold_for_plots, 'ASD-gs', isLog);
% fits.ASD_gs.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
%         data.nt, fold_for_plots, fits.ASD_gs.lbl, sc);
% fig_svfcn(fits.ASD_gs.fig, 'decision-ASD-gs', 'png');

[fits.ML, wf, sc] = reg.cvMaxScoreGrid(data, [nan nan nan], mlFcn, {}, ...
    M.rsqFcn, {}, foldinds, 'ML', 1);
fits.ML.fig = plot.prepAndPlotKernel(data.Xxy, wf, data.ns, ...
    data.nt, fold_for_plots, fits.ML.lbl, sc);
fig_svfcn(fits.ML.fig, 'decision-ML', 'png');

fits.isLinReg = isLinReg;
io.updateStruct(dat_fnfcn('decision'), fits);
