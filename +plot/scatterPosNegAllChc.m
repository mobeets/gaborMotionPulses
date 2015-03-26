function scatterPosNegAllChc(v)
    % next, only display responses to small stimstrength
    
    inds = v.dirprob == 0;
    if sum(inds) < 100
        disp('using dirstrength');        
        inds = v.dirstrength < 10;
        disp(['mean dirstrength= ' num2str(mean(v.dirstrength(inds)))]);
    end
    disp(['showing ' num2str(sum(inds)) ' trials']);
    disp(['ignoring ' num2str(sum(~inds)) ' trials']);

    hold on;
    chc = v.C;
    A = chc&inds;
    B = ~chc&inds;
    e1 = randn(size(v.Y(A)))*.1;
    e2 = randn(size(v.Y(B)))*.1;
%     scatter3(v.Ypos(A), v.Yneg(A), v.Yh(A)+e1, 50, [1 0.5 0.5]);
%     scatter3(v.Ypos(B), v.Yneg(B), v.Yh(B)+e2, 50, [0.5 0.5 1]);
    scatter3(v.Ypos(A), v.Yneg(A), v.Y(A)+e1, 50, [1 0.5 0.5], 'filled');
    scatter3(v.Ypos(B), v.Yneg(B), v.Y(B)+e2, 50, [0.5 0.5 1], 'filled');
    xlabel('pos'), ylabel('neg'), zlabel('all');
    title([v.name '   Yposneg-corr=' num2str(v.corr_Yposneg) ...
        ' sc=' num2str(v.score) '  |wpos|=' num2str(v.wpos_norm) ...
        '  |wneg|=' num2str(v.wneg_norm)])
end
