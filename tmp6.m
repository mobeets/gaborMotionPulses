dt = '20150304a';
d = io.loadDataByDate(dt, true);
fs = io.loadFitsByDate(dt, 'data/evirepb-nancy/fits');
% fs = io.loadFitsByDate(dt, 'data/ml-nancy/fits');
% f0 = fs.decision;
f0 = fs.LIP_1;

foldinds = f0.foldinds{end};
isLinReg = f0.isLinReg{end};
X = d.X;
D = d.D;

f = f0.ASD{end};
% f = f0.ML{end};

% Y = d.R;
Y = d.Y_all(:,1);


gs = io.loadFitsByDate(dt, 'data/refactor-nancy/fits');

%%
y = load('tmp.mat');
foldinds = y.foldinds;

%%
clear obj
obj.foldinds = foldinds;
% scoreObj = reg.getScoreObj(isLinReg, 'pctCorrect');
scoreObj = reg.getScoreObj(isLinReg, 'rsq');
% obj = reg.getObj_ML(X, Y, obj);
obj = reg.getObj_ASD(X, Y, D, scoreObj, obj);
obj = reg.fitAndScore(X, Y, obj, scoreObj);
% obj = tools.rmfieldsRegexp(obj, {'Fcn$', 'FcnArgs$'}, true);


