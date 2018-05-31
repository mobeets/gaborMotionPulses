function obj = getSeparableRF(w)
    
    [U, s, V] = svd(w, 0);
    S = diag(s);
    separabilities = S/norm(S,1);
    
    clear obj;
    obj.separability_index = separabilities(1); % norm-1
    obj.inseparability_index = 1 - S(1)/norm(S,2); % norm-2 via Linden
    
    wf_svdf = @(jj) U(:,jj)*S(jj)*V(:,jj)';
    obj.wfSvd_U = U;
    obj.wfSvd_V = V;
    obj.wfSvd_S = S;
    obj.wfSvd_1 = wf_svdf(1);
    obj.wfSvd_2 = wf_svdf(2);
    obj.wfSvd_3 = wf_svdf(3);
    
    obj.spatial_RF = U(:,1);    
    obj.temporal_RF = V(:,1);
    obj.scalar_RF = S(1);
    
end
