function plotAndSaveKernel(obj, data, figdir, add_dt, add_score)
% 
% plots the kernel and saves to png
% 
    if nargin < 5
        add_score = false;
    end
    if nargin < 4
        add_dt = false;
    end
    fig_lblfcn = @(lbl, sc) [lbl ' sc=' num2str(sprintf('%.2f', sc))];
    
    if ~isempty(figdir)
        wf = obj.mu;
        if prod(obj.shape) < numel(wf)
            wf = wf(1:end-1);
        end
        tags = strsplit(obj.label, '-');
        isCell = numel(tags) > 2;
        t1 = median(data.stim.targ1XY);
        t2 = median(data.stim.targ2XY);
        if isCell && data.neurons{str2num(tags{2})}.targPref == 2
            wf = -wf;
            t3 = t1;
            t1 = t2;
            t2 = t3;
        end
        lbl = obj.label;
        if add_dt
            lbl = [obj.dt '-' obj.label];
        end
        if add_score
            lblS = [lbl ', score = ' sprintf('%.2f', obj.score)];
        else
            lblS = lbl;
        end
        
%         if ~isCell
%             t1 = nan;
%             t2 = nan;
%         end

        plot.plotKernel2(reshape(wf, obj.shape(1), obj.shape(2)), ...
            data.Xxy, nan, [0 0], t1, t2, lblS); 
        plot.saveFigure(lbl, figdir);
    end
end
