%% load fit data

% fitdir = '20150615';
% fitdir = '20150716';
fitdir = 'new';
fitbasedir = 'data';

% ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
% vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
% vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
% vu = [vp vn];
fitdir = fullfile('data', 'fits', fitdir);
vu = tools.makeFitSummaries(fitdir, true, 'ASD');
vu = figure.filterData(vu);

% mark RFs that are not uniform in space
tm = arrayfun(@(v) var(v.wfSvd_1(:,1)), vu, 'uni', 0);
[vu.wfVar] = tm{:};
vuMT = vu([vu.isMT]);
vuMT_var = vuMT(log([vuMT.wfVar])>-17);
% vuMT_var = vuMT;
%%

figBaseDir = '~/code/gaborMotionPulses/data/figs/';
dirNames = arrayfun(@(x) dir(fullfile(figBaseDir, ...
    ['Figure0' num2str(x) '*'])), 1:5, 'uni', 0);
figDirFcn = @(num) fullfile(figBaseDir, dirNames{num}.name);

%%

figDir = '~/code/gaborMotionPulses/data/figs/';
figure.figure2;

%%

% methods
% figure1

% ASD methods with example fits
figDir = figDirFcn(2);
figure.figure2;

% ASD/hyperflow overlay example cells
figDir = figDirFcn(3);
figure.figure3;

% pairwise noise correlation and CP vs. RF
figDir = figDirFcn(4);
figure.figure4;

% decoding
figure5_reGenData = false;
figure.figure5;
