%% load fit data

fitdir = '20150615';
fitbasedir = 'data';
ff = @(mnkNm) fullfile(fitbasedir, [fitdir '-' mnkNm], 'fits');
vp = tools.makeFitSummaries(ff('pat'), false, 'ASD');
vn = tools.makeFitSummaries(ff('nancy'), true, 'ASD');
vu = [vp vn];

%% filter fit data

nMT = sum([vu.isMT]);

% ignore sessions where monkey's pctCorrect < 75%
vu = vu([vu.pctCorrect] >= 0.75);

vuMT = vu([vu.isMT]); nMT = numel(vuMT);
% ignore cells with dPrime < 0.4
vuMT = vuMT([vuMT.dPrime] >= 0.4);
% ignore cells with ntrials < 100
vuMT = vuMT([vuMT.ntrials] >= 100);

warning(['Removing ' num2str(nMT-numel(vuMT)) ' of ' ...
    num2str(nMT) ' MT cells.']);

%%

figBaseDir = '~/Dropbox/gaborMotionPulseASD/figures/';
dirNames = arrayfun(@(x) dir(fullfile(figBaseDir, ...
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
