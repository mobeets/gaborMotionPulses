nx = 3;
dx = 1:nx;
D = asd.sqdist.time(nx, 1);

w = 5*rand(1,nx)';
N = 10;
ws = ones(N,1)*w';
ws(:,2) = linspace(0.1, 3, N);
ds = nan(N,1);
ros = nan(N,1);

clrs = summer(N);
figure(1); clf;

for ii = 1:N
    w = ws(ii,:)';
    
    % d | w, ro, D
    
    ro = 0.2;
    pr = @(x) asd.prior(ro, D, x);

    P = @(x) mvnpdf(w, 0*dx', pr(x));
    d = fmincon(@(x) -log(P(x)), 1, [],[],[],[], 0.4, 100);
    ds(ii) = d;    

    xs = linspace(0.4, d+5, 1e3);
    ys = arrayfun(P, xs);
%     ys = ys./sum(ys);
    ys = ys./max(ys);
    
    % using expected value and mode to estimate shape, scale for gamma
    mu = xs*ys'/sum(ys);
    mo = d;
    k = mu/(mu-mo); th = mu/k; % real calculation, but this has error
    ys0 = gampdf(xs, k, th);
    ys0 = ys0./max(ys0);
%     ys0 = ys0./sum(ys0);
    
    subplot(2,1,1); hold on;    
    plot(xs, ys, '--', 'color', clrs(ii,:));
    plot(xs, ys0, 'color', clrs(ii,:), 'HandleVisibility','off');
    xlabel(['delta | ro = ' num2str(ro)]);
    
    % ro | w, d, D

    d = 0.5;
    pr = @(ro) asd.prior(ro, D, d);

    P = @(x) mvnpdf(w, 0*dx', pr(x));
    ro = fmincon(@(x) -log(P(x)), 1, [],[],[],[], -20, 20);
    ros(ii) = ro;

    xs = linspace(-20, 20, 1e3);
    ys = arrayfun(P, xs);
%     ys = ys./sum(ys);
    ys = ys./max(ys);
    xs = exp(xs);
    
    % using expected value and mode to estimate shape, scale for gamma
    mu = xs*ys'/sum(ys);
    mo = d;
    k = mu/(mu-mo); th = mu/k; % real calculation, but this has error
    ys0 = gampdf(xs, k, th);
    ys0 = ys0./max(ys0);
%     ys0 = ys0./sum(ys0);
    
    subplot(2,1,2); hold on;    
    plot(xs, ys, '--', 'color', clrs(ii,:));
    plot(xs, ys0, 'color', clrs(ii,:), 'HandleVisibility','off');
    xlabel(['exp(ro) | delta = ' num2str(d)]);
end

subplot(2,1,1);
legend(arrayfun(@(x) sprintf('%0.2f', x{1}), num2cell(ws(:,2)), 'uni', 0));

[ros ds ws(:,1) ws(:,2)]

% xs = sqrt(sum(abs(ws).^2,2));
xs = ws(:,2);
figure(2); hold on;
subplot(2,1,1); hold on;
scatter(xs, ds, 'filled');
ylabel('mode of delta | w_1, w_2, ro');
subplot(2,1,2); hold on;
scatter(xs, exp(ros), 'filled');
xlabel(['w_1 = ' num2str(w(1)) ', w_2 = x-axis']);
ylabel('mode of exp(ro) | w_1, w_2, delta');
