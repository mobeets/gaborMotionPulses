function cell = RfCenterOfMass(cell)

    % normalized stimulus location
    xs = zscore(cell.Xxy);

    % thresholded spatial RF
    if isfield(cell, 'wsep')
        ys0 = cell.wsep.spatial_RF;
    else
        warning('Cannot use wfSvd_1 for centerOfMass. Using w instead.');
        ys0 = cell.w;
    end
    ys = ys0*sign(sum(ys0)); % make positive
    ys(ys < 0) = 0; % ignore subfields
    ys = ys / sum(ys); % normalize so we can use as weights
    
    % mean of location weighted by RF strength
    xc = (xs(:,1)'*ys);
    yc = (xs(:,2)'*ys);
    
    [theta, rho] = cart2pol(xc, yc);
    cell.rf_center = [xc yc];
    cell.rf_ecc = rho;
    cell.rf_theta = theta;
    
    warning('Need to validate new center of mass estimate...');
    ys = ys0;
    if (cell.targPref ~= 1) % used to be == 1 but that gave problems
        % flip weights so that < 0 always refers to subfield
        ys = -ys;
    end
    ps = cell.Xxy; % spatial positions
    ps = zscore(ps); % z-scored positions so we can compare across sessions
    isSubfield = ys < 0; % ignore mass and positions of subfield
    ps = ps(~isSubfield,:);
    ys = ys(~isSubfield);
    pts = mean(bsxfun(@times, ps, ys));
    xc = pts(1); yc = pts(2);
    [theta, rho] = cart2pol(xc, yc);
    cell.rf_center = [xc yc];
    cell.rf_ecc = rho;
    cell.rf_theta = theta;
end
