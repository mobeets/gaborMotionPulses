function [ps, nms] = spikeRateByTrial(dt, rsq_thresh)
    if nargin < 1
        dts = {'20130502', '20130514', '20130515', '20130517', ...
            '20130611', '20140213', '20140218', '20140226', '20140303', ...
            '20140304', '20140305', '20140306', '20140307', '20140310'};
    else
        dts = {dt};
    end
    if nargin < 2
        rsq_thresh = 0.1;
    end
    ps = [];
    nms = [];
    for ii = 1:numel(dts)
        dt = dts{ii};
        d = io.loadDataByDate(dt);
        ncells = size(d.Y_all, 2);
        for jj = 1:ncells        
            ys = d.Y_all(:,jj)';
            ys = ys(~isnan(ys));
            xs = 1:numel(ys)';
            p = polyfit(xs, ys, 1);
            yh = polyval(p, xs);
            rsq = 1 - sum((ys-yh).^2)/((numel(ys)-1)*var(ys));
            nm = [dt '-' num2str(jj)];
            
            ps = [ps; [p rsq]];
            nms = [nms {nm}];

            if rsq > rsq_thresh
                figure; hold on;                
                disp(nm);
                lbl = [nm ' m=' num2str(p(1)) ' rsq=' num2str(rsq)];
                title(lbl);
                plot(xs, yh);
                plot(smooth(ys));
            end
        end
    end
end
