%% load fit data

fitdir = '20150615';
fitbasedir = 'data';
ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
vu = [vp vn];
vuMT = vu([vu.isMT]);

%%

figBaseDir = '~/Dropbox/gaborMotionPulseASD/figures/';
dirNames = arrayfun(@(x) dir(fullfile(fdir, ...
    ['Figure0' num2str(x) '*'])), 1:5, 'uni', 0);
figDirFcn = @(num) fullfile(figBaseDir, dirNames{num}.name);

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
