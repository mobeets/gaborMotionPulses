dt = '20150304a';
d = io.loadDataByDate(dt, true);
fs = io.loadFitsByDate(dt, 'data/evirepb-nancy/fits');
% fs = io.loadFitsByDate(dt, 'data/ml-nancy/fits');
f0 = fs.decision;
% f0 = fs.MT_5;

foldinds = f0.foldinds{end};
isLinReg = f0.isLinReg{end};
X = d.X;
D = d.D;

f = f0.ASD{end};
% f = f0.ML{end};

Y = d.R;
% Y = d.Y_all(:,5);

%%
y = load('tmp.mat');
foldinds = y.foldinds;

%%
clear obj
obj.foldinds = foldinds;
scoreObj = reg.getScoreObj('pctCorrect', isLinReg);
% scoreObj = reg.getScoreObj('rsq', isLinReg);
% obj = reg.getObj_ML(X, Y, obj);
obj = reg.getObj_ASD(X, Y, D, scoreObj, obj);
obj = reg.fitAndScore(X, Y, obj, scoreObj);
obj = tools.rmfieldsRegexp(obj, {'Fcn$', 'FcnArgs$'}, true);
