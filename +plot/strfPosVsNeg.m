function vals = strfPosVsNeg(fitstr)
    if nargin < 1
        fitstr = 'ASD';
    end
    dts = io.getDates('fits');
    c = 0;
    for ii = 1:numel(dts)
        dt = dts{ii};
        d = io.loadDataByDate(dt);
        fs = io.loadFitsByDate(dt);
        nms = fieldnames(fs);
        for jj = 1:numel(nms)
            fit = fs.(nms{jj}).(fitstr){end};            
            wf = fit.mu(1:end-1);
            mu = reshape(wf, d.ns, d.nt);
            
            muPos = nan(size(mu));
            muPos(mu > 0) = mu(mu > 0);
            muNeg = nan(size(mu));
            muNeg(mu < 0) = mu(mu < 0);

            c = c+1;
            
            vals(c).type = nms{jj}(1:2);
            
            vals(c).wf = wf;            
            vals(c).wf_mean = mean(wf);
            vals(c).wf_sdev = std(wf);
            
            vals(c).mu = mu;
            vals(c).mu_mean = mean(mu);
            vals(c).mu_sdev = std(mu);
            
            vals(c).muPos = muPos;
            vals(c).muPos_mean = nanmean(muPos);
            vals(c).muPos_sdev = nanstd(muPos);
            
            vals(c).muNeg = muNeg;
            vals(c).muNeg_mean = nanmean(muNeg);
            vals(c).muNeg_sdev = nanstd(muNeg);
        end
    end
end
