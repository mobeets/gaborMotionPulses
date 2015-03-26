% 
% 
%
%
%
%
%

%%
rng(1234);
ntrials = 1000;
C = randn(ntrials,1) > 0;

gA = 1;
gB = -0.1;
vA = 5;
vB = 0.5;

eA = sqrt(vA)*randn(ntrials,1);
eB = sqrt(vB)*randn(ntrials,1);
eZ = sqrt(vA+vB)*randn(ntrials,1);
A = gA*C + eA;
B = gB*C + eB;
% Z = (gA+gB)*C + eZ;
Z = (gA+gB)*C + eA + eB;



CPA = tools.AUC(A(C), A(~C));
CPB = tools.AUC(B(C), B(~C));
CPZ = tools.AUC(Z(C), Z(~C));
disp('---');
disp('     ')
cPA = corrcoef(A, +C);
cPB = corrcoef(B, +C);
cPZ = corrcoef(Z, +C);
disp('corrZ            corrA          corrB');
disp(num2str([cPZ(2) cPA(2) cPB(2)]));
disp('     ')
disp('Z             A           B');
disp(num2str([CPZ CPA CPB]));

%%
rng(1234);
ntrials = 1000;
C = randn(ntrials,1) > 0;

Ypos = C + 0.5*randn(ntrials,1);
Yneg = C + 1*randn(ntrials,1);
Yres = C + 1*randn(ntrials,1);
Y = Ypos+Yneg+Yres;

CP_pos = tools.AUC(Ypos(C), Ypos(~C));
CP_neg = tools.AUC(Yneg(C), Yneg(~C));
CP_res = tools.AUC(Yres(C), Yres(~C));
CP_all = tools.AUC(Y(C), Y(~C));
disp('---');
disp('CP_all      CP_pos      CP_neg      CP_res');
disp(num2str([CP_all CP_pos CP_neg CP_res]));

%%

rng(1234);

ntrials = 1000;
ngabors = 2;
% X = repmat(randn(ntrials, 1), 1, 2);
Yneg = randn(ntrials, ngabors);

e = randn(ntrials,1);
Yres = Yneg*w + e;

thresh = 0;
eps = 1.5*randn(ntrials,1);
% figure; hist(Y);
C = (Yres + eps) > thresh;

wpos = w;
wpos(wpos < 0) = 0;
wneg = w;
wneg(wneg > 0) = 0;

Ypos = Yneg*wpos;
Yneg = Yneg*wneg;

% close all;
figure(3); hold all
scatter(Ypos, Yneg);

% scatter(X(:,1), X(:,2));

% scatter(X(:,2), C);

disp('-----');

CP_Y = tools.AUC(Yres(C), Yres(~C));

Y_r = Yres - Ypos - Yneg;
CP_Y_r = tools.AUC(Y_r(C), Y_r(~C));

Y_rneg = Yres - Yneg;
CP_Y_rneg = tools.AUC(Y_rneg(C), Y_rneg(~C));

Y_rpos = Yres - Ypos;
CP_Y_rpos = tools.AUC(Y_rpos(C), Y_rpos(~C));

disp(['w1           w2  CP_all  CP_rneg     CP_rpos     CP_r']);
disp(num2str([w', CP_Y, CP_Y_rneg, CP_Y_rpos, CP_Y_r]));
