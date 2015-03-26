function wrs = wposnegRatio(wpos, wneg)
    wrs = nan(numel(wpos),1);
    for ii = 1:numel(wpos)
        top = norm(wneg{ii},1);
        bot = norm(wpos{ii},1);
        wrs(ii) = top/bot;
    end
end
