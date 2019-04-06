function [bp, mu, sigma] = gauss2dcirc(data, sigMult, sigma)
% data [n x 2]
% plot(bp(1,:), bp(2,:)) will plot a circle
% 
    if nargin < 2
        sigMult = 1;
    end
    if nargin < 3
        assert(size(data,2)==2);
        mu = nanmean(data)';
        sigma = nancov(data);
    else
        mu = zeros(size(sigma,1),1);
    end
    
    nd = size(sigma,1);
    if nd == 2
        tt = linspace(0, 2*pi, 60)';
        x = cos(tt); y = sin(tt);
        ap = [x(:) y(:)]';
    else
        % https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates
        ts = repmat({linspace(0, pi, 30)}, nd-1, 1);
        ts{end+1} = linspace(0, 2*pi, 30);
        T = cell(1, numel(ts));
        [T{:}] = ndgrid(ts{:});
        phis = nan(nd, numel(T{1}(:)));
        for ii = 1:nd
            phis(ii,:) = T{ii}(:);
        end
        ap = nan(size(phis));
        for ii = 1:nd
            if ii == 1
                ap(1,:) = cos(phis(1,:));
            elseif ii < nd
                ap(ii,:) = prod(sin(phis(1:ii-1,:))).*cos(phis(ii,:));
            else
                ap(ii,:) = prod(sin(phis(1:ii,:)));
            end
        end
    end
    [v,d] = eig(sigma);
    d(abs(d) < eps) = 0; % sometimes we get slightly negative eigs
    d = sigMult*sqrt(d);
    bp = (v*d*ap) + repmat(mu, 1, size(ap,2)); 

end
