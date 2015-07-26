function cleanupForPrint(varargin)
% clean up figure for print
% cleanupForPrint(figure_handle, varargin)
% Argument Options
%   'FontSize'     [10] (default)
%   'PaperSize'    [174 100] (units: mm)
% 
if nargin < 1 || ~isa(varargin{1}, 'matlab.ui.Figure')
    h = gcf;
    if isempty(varargin)
        args = {};
    else
        args = varargin(:);
    end
else
    h = varargin{1};
    if numel(varargin)==1
        args = {};
    else
        args = varargin(2:end);
    end
end

p = inputParser();
p.addOptional('PaperSize', [174 100])
p.addOptional('PaperUnits', 'pixels')
p.addOptional('FontSize', 10)
p.addOptional('FontName', 'Helvetica')
p.addOptional('FontWeight', 'normal')
p.addOptional('Linewidth', .8)
p.addOptional('TickDir', 'out')
p.addOptional('Box', 'off')
p.addOptional('plot2svg', false);
p.parse(args{:})

opts = p.Results;
opts.PaperSize     = opts.PaperSize/10;
opts.PaperPosition = [0*opts.PaperSize opts.PaperSize];
options = fieldnames(opts);

set(0, 'currentfigure', h)

fc = get(h,'children');

% Loop through all the children
for kChild = 1:length(fc)
    
    switch class(fc(kChild))
        case 'matlab.graphics.illustration.Legend'
            set(fc(kChild), 'Box', 'off')
        
        case 'matlab.graphics.axis.Axes'
            
            set(gcf,'currentaxes',fc(kChild))
            xt = get(fc(kChild), 'Xtick');
            xlim(xt([1 end]));
            yt = get(fc(kChild), 'Ytick');
            ylim(yt([1 end]));
            
            
            ac = get(fc(kChild),'children');
            ap = get(fc(kChild));
            
            for ii = 1:numel(options)
                if isfield(ap, options{ii})
                    set(fc(kChild), options{ii}, opts.(options{ii}));
                end
            end
            
            % And loop through all the children of each axis:
            for kAxis = 1:length(ac)
                af = get(ac(kAxis));
                
                % Loop through axis options and modify the axis
                for ii = 1:numel(options)
                    if isfield(af, options{ii}) && ~strcmpi(options{ii}, 'Linewidth')
                        set(ac(kAxis), options{ii}, opts.(options{ii}));
                    end
                end
            end
            
            ht = get(fc(kChild),'title');
            set(ht,'FontName',opts.FontName,'FontSize',opts.FontSize);
            hx = get(fc(kChild),'xlabel');
            set(hx,'FontName',opts.FontName, 'FontWeight', opts.FontWeight,'FontSize',opts.FontSize);
            hy = get(fc(kChild),'ylabel');
            set(hy,'FontName',opts.FontName, 'FontWeight', opts.FontWeight,'FontSize',opts.FontSize);
            %         hl  = findobj(gcf,'Type','axes','Tag','legend');
            %         set(hl,'box','off')
    end
    
end


set(h, 'PaperSize', opts.PaperSize, 'PaperPosition', opts.PaperPosition)
set(h, 'Color', 'w')

% pixel hack for plot2svg
if opts.plot2svg
    set(h, 'Units', 'centimeters');
    pos = get(h, 'Position');
    pos(3) = opts.PaperSize(1); % Select the width of the figure in [cm]
    pos(4) = opts.PaperSize(2); % Select the height of the figure in [cm]
    set(h, 'Position', pos);
    set(h, 'PaperType', 'a4letter');
end