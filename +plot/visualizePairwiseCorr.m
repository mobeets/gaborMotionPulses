function S = visualizePairwiseCorr(r, c)
% visualize joint responses for neuron pairs
% S = vizualizePairwiseCorr(spikeCount, Condition)
% INPUT:
%   spikeCount = [nTrials x 2]
%   Condition  = [nTrials x 1] logical condition

cmap = lines;

ix0=~any(isnan(r),2);
ix1=~any(isnan(c),2);
ix = ix0&ix1;
r = r(ix,:);
c = c(ix);

% get group mean and covariances
mu1 = mean(r(c, :)); sig1 = cov(r(c, :));
mu2 = mean(r(~c, :)); sig2 = cov(r(~c, :));
mu  = mean(r);

scatter(r(c,1), r(c,2), 5, cmap(2,:)); hold on
scatter(r(~c,1), r(~c,2), 5, cmap(1,:));

plotellipse(mu1, sig1, 1, 'Color', cmap(2,:));
plotellipse(mu2, sig2, 1, 'Color', cmap(1,:));

% compute Fisher Linear Discriminant
fld = (sig1+sig2)\(mu2-mu1)'; % standard (use full covariance matrix)
fld0 = (diag(diag(sig1))+diag(diag(sig2)))\(mu1-mu2)'; % correlation blind (zero off-diagonal values)

% line functions
q   = @(x) (fld(2)/fld(1))*x + mu(2) - (fld(2)/fld(1))*(mu(1)+sig1(1)/3);
qOr = @(x) (-fld(1)/fld(2))*x + (mu(2) - (-fld(1)/fld(2))*mu(1));
q1 = @(x) (-fld0(1)/fld0(2))*x + (mu(2) - (-fld0(1)/fld0(2))*mu(1));

% find intersection of FLD and FLDorth (why did I do it this way?)
% Jay, if you get bored... simplify this shit!
a0 = (fld(2)/fld(1));
b0 = -1/a0;
c0 = mu(2) - (fld(2)/fld(1))*(mu(1)+sig1(1)/3);
d0 = (mu(2) - (-fld(1)/fld(2))*mu(1));
x0 = (d0-c0)/(a0-b0);

% find optimal criterion for the 2 cases
crit1 = (mu1+mu2)*fld/2;
crit0 = (mu1+mu2)*fld0/2;

% project r and decode under the two conditions
rproj = r*fld;
S.mu1 = mu1;
S.mu2 = mu2;
S.sig1 = sig1;
S.sig2 = sig2;
S.fld  = fld;
S.fldCorrBlind = fld0;
S.pcOpt     = 1-mean(c == (rproj > crit1));
S.pcBlind   = mean(c == (r*fld0 > crit0));
S.critOpt   = crit1;
S.critBlind = crit0;




qd = ezplot(q, [0 max(xlim) 0 max(ylim)]);
qo = ezplot(qOr, [0 max(xlim) 0 max(ylim)]);
set(qo, 'Linewidth', 2, 'Color', 'k', 'Linestyle', '-');
set(qd, 'Linewidth', 2, 'Color', 'k', 'Linestyle', '-');
qe = ezplot(q1, [0 max(xlim) 0 max(ylim)]);
set(qe, 'Linewidth', 2, 'Color', 'k', 'Linestyle', '--');

% plot the projection
rproj = zscore(rproj-crit1); % zscore to scale
sc = mean(std(r)); % scale value
rproj = rproj*sc; %scale projection
bins = linspace(min(rproj), max(rproj), 25); % binning
cnt1 = histc(rproj(c), bins)'; 
cnt2 = histc(rproj(~c), bins)';
cnt1 = cnt1/sum(cnt1)*sc*10; % scale histogram
cnt2 = cnt2/sum(cnt2)*sc*10; % scale histogram
% stretch to have jagged edges to look like a bar plot
bins = circshift(reshape([bins; bins], [], 1),1);
cnt1 = reshape([cnt1; cnt1], [], 1); cnt1([1 2 end]) = 0;
cnt2 = reshape([cnt2; cnt2], [], 1); cnt2([1 2 end]) = 0;
% build rotation matrix
theta = cart2pol(-fld(1), fld(2));
s = sign(theta);
R = [-cos(theta) sin(theta); sin(theta) cos(theta)];
X1 = [bins(:) s*cnt1(:)]*R;
X2 = [bins(:) s*cnt2(:)]*R;

ab = x0;
plot(X1(:,1)+ab, q(ab)+ X1(:,2), 'Color', cmap(2,:));
plot(X2(:,1)+ab, q(ab)+ X2(:,2), 'Color', cmap(1,:));
xlabel('Neuron 1')
xlabel('Neuron 2')
title('')

function [h,el] = plotellipse(mu, covmat, r, varargin)
% h = plotellipse(mu, covmat, stdev, varargin);
% plots an ellipse of constant probability (i.e. a contour of
% constant deviation stdev) for a given bivariate gaussian distr.
% Inputs:
%    mu = column vector with the mean of the distr.
%    covmat = 2x2 covariance matrix
%    stdev = the standard deviation contour we wish to plot (e.g. 1, 2, .2, etc)
% Output: 
%    h = handle to plotted line
%    el = 100x2 matrix where the two columns provide the x and y values of the 
%         ellipse

[U, S, V] = svd(covmat);
thet = [0:(2*pi)/99:(2*pi+.0001)];
el = [r*U*sqrt(S)*[cos(thet); sin(thet)]]';
el(:,1) = el(:,1)+mu(1);
el(:,2) = el(:,2)+mu(2);
h = plot(el(:,1), el(:,2), varargin{:});