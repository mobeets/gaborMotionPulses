function vals = pairwiseCorr(vs)
    vals = struct([]);

    dts = unique({vs.dt});
    for ii = 1:numel(dts)
        dt = dts{ii};
        inds = strcmp(dt, {vs.dt});
        vsc = vs(inds);
        Y = [vsc.Y];
        Yres = [vsc.Yres];
        nms = {vsc.dt};
        tps = {vsc.type};
        scs = [vsc.score];
        if size(Y,2) < 2
            continue;
        end
        pairs = nchoosek(1:size(Y,2), 2);
        for jj = 1:size(pairs,1)        
            a = pairs(jj,1);
            b = pairs(jj,2);
            
            Y1 = Y(:,a);
            Y2 = Y(:,b);
            ix = ~isnan(Y1) & ~isnan(Y2);
            p0 = corr(Y1(ix), Y2(ix));
            
            Y1 = Yres(:,a);
            Y2 = Yres(:,b);
            ix = ~isnan(Y1) & ~isnan(Y2);
            p1 = corr(Y1(ix), Y2(ix));
            
            p2 = corrcoef(vsc(a).mu, vsc(b).mu);        

            val.nm1 = nms{a};
            val.nm2 = nms{b};
            val.spikeCorr = p0;
            val.spikeResidCorr = p1;
            val.rfCorr = p2(2);
            if strcmp(tps(a), tps(b))
                if strcmp(tps(a), 'MT')
                    val.type = 'MT-MT';
                else
                    val.type = 'LIP-LIP';
                end
            else
                val.type = 'MT-LIP';
            end
            vals = [vals val];
        end
    end
end
