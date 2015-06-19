function plotAndSaveKernel(obj, data, figdir, add_dt)
% 
% plots the kernel and saves to png
% 
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
        if isCell && data.neurons{str2num(tags{2})}.targPref == 2
            wf = -wf;
        end
        if add_dt
            lbl = [obj.dt '-' obj.label];
        else
            lbl = obj.label;
        end
        t1 = median(data.stim.targ1XY);
        t2 = median(data.stim.targ2XY);
        if ~isCell
            t1 = nan;
            t2 = nan;
        end

        plot.plotKernel2(reshape(wf, obj.shape(1), obj.shape(2)), ...
            data.Xxy, nan, [0 0], t1, t2, lbl);
%         [fig, ha] = plot.plotKernel(data.Xxy, ...
%             reshape(wf, data.ns, data.nt), ...
%             nan, fig_lblfcn([obj.dt '-' obj.label], obj.score));
%         
%         axes(ha(1)); hold on;
%         scatter(0, 0, 50, [0.8 0.2 0.2], 'filled');
%         t1 = median(data.stim.targ1XY);
%         t2 = median(data.stim.targ2XY);
%         scatter(t1(1), t1(2), 50, [0.2 0.8 0.2], 'filled');
%         scatter(t2(1), t2(2), 50, [0.2 0.6 0.2], 'filled');
%         delta = 0.1;
%         vs = [t1; t2; 0 0];
%         vs = [vs-delta; vs+delta];
%         vs = [vs; data.stim.gaborXY];
%         vsb = [min(vs); max(vs)];
%         lb = min(vsb(:));
%         ub = max(vsb(:));
%         for ii = 1:numel(ha)
%             axes(ha(ii));
%             xlim([lb ub]);
%             ylim([lb ub]);
%         end
%         
        plot.saveFigure(lbl, figdir);
    end
end
