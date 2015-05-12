function fig = quickPmfByDate(dt, isNancy, nbins, byPulse)
    if nargin < 4
        byPulse = false;
    end
    if nargin < 3
        nbins = 10;
    end
    data = io.loadDataByDate(dt, isNancy);
    Y = data.R;
        
    fig = figure; hold on; set(gcf,'color','w');
    xlabel('marginal stimulus strength');
    ylabel('% pref choice');
    title(dt);
    
    X = sum(data.X, 2);
    mx = max(abs(X));
    dv = floor(mx/(nbins/2));
    M = [floor(X/dv) Y];
    plotPctByGroup(M, 'k');
    
    % pmf for each temporal pulse
    if byPulse
        X = reshape(data.X, numel(Y), data.ns, data.nt);
        clrs = gray(size(X,3)+1);
        for ii = 1:size(X,3)
            Xt = sum(X(:,:,ii),2);
            mx = max(abs(Xt));
            dv = floor(mx/(nbins/2));
            M = [floor(Xt/dv) Y];
            plotPctByGroup(M, clrs(ii,:));
        end
    end
end

function plotPctByGroup(M, clr)
    [ii,~,jj]  = unique(M(:,1), 'rows');
    ps = accumarray(jj, M(:,2), [], @mean);
    ns = accumarray(jj, M(:,2), [], @numel);
    errs = sqrt(ps.*(1-ps)./ns);
    errorbar(ii, ps, errs, 'r', 'LineStyle', 'none', 'LineWidth', 2);
    scatter(ii, ps, 'ko', 'MarkerFaceColor', clr);
    plot(ii, ps, 'Color', clr, 'LineWidth', 2);
end
