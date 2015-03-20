function vals = strfPosVsNeg(vals)
% 
%   vals - from plot.goodFits()
% 
    inds = [vals.score]./[vals.scoreSdev] > 1;
    vals = vals(inds);
    for ii = 1:numel(vals)
        mu = vals(ii).mu;
        muPos = nan(size(mu));
        muPos(mu > 0) = mu(mu > 0);
        muNeg = nan(size(mu));
        muNeg(mu < 0) = mu(mu < 0);

        wf = mu(:);
        vals(ii).wf = wf;            
        vals(ii).wf_mean = mean(wf);
        vals(ii).wf_sdev = std(wf);

        vals(ii).mu = mu;
        vals(ii).mu_mean = mean(mu);
        vals(ii).mu_sdev = std(mu);

        vals(ii).muPos = muPos;
        vals(ii).muPos_sum = nansum(muPos);
        vals(ii).muPos_mean = nanmean(muPos);
        vals(ii).muPos_sdev = nanstd(muPos);

        vals(ii).muNeg = muNeg;
        vals(ii).muNeg_sum = nansum(muNeg);
        vals(ii).muNeg_mean = nanmean(muNeg);
        vals(ii).muNeg_sdev = nanstd(muNeg);
    end
end
