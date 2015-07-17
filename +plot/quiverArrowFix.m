function quiverArrowFix(qh, headWidth, headLength, varargin)
% quiverArrowFix(qh)
%
% replaces quiver arrows with annotation arrows
%
% qh - quiver figure handle
%     [x,y] = meshgrid(0:0.2:2,0:0.2:2);
%     u = cos(x).*y; v = sin(x).*y;
%     plot.quiverArrowFix(quiver(x,y,u,v), 'LineWidth', 2);
%  OR:
%     plot.quiverArrowFix(quiver(x,y,u,v), 30, 80, 'LineWidth', 2);
%
% src: http://stackoverflow.com/questions/18776172/in-matlab-how-do-i-change-the-arrow-head-style-in-quiver-plot
%

if ischar(headWidth)
    if nargin < 4
        varargin = {headWidth, headLength};
    else
        varargin = [{headWidth, headLength} varargin];
    end
end
if nargin < 3 || ischar(headWidth)
    headLength = 30;
end
if nargin < 2 || ischar(headWidth)
    headWidth = 80;
end

set(qh, 'Visible', 'off');
% correct for different graphics engine in newer matlabs
if datenum(version('-date')) >= datenum('September 15, 2014')
    x  = qh.XData;
    y  = qh.YData;
    dx = qh.UData;
    dy = qh.VData;
    
    arwlngth = max(sqrt(dx.^2 + dy.^2));
    dx = 1.5*dx/arwlngth;
    dy = 1.5*dy/arwlngth;
%     rnge = max(x)-min(x);
    
    clr = [0 0 0];
    pos = get(gca, 'position');
    xd = xlim;
    yd = ylim;
    
    x1 = x - dx;
    x2 = x + dx;
    
    y1 = y - dy;
    y2 = y + dy;
    
    % coordinate transform 1
    y1 = (y1-yd(1))/diff(yd);
    y2 = (y2-yd(1))/diff(yd);
    x1 = (x1-xd(1))/diff(xd);
    x2 = (x2-xd(1))/diff(xd);
    % rescale to gca
    x1 = x1*pos(3) + pos(1);
    x2 = x2*pos(3) + pos(1);
    y1 = y1*pos(4) + pos(2);
    y2 = y2*pos(4) + pos(2);
    
    
    
    for ii = 1:numel(x1)
        % ah = annotation('arrow',[x1(ii) x2(ii)], [y1(ii) y2(ii)], ...
        ah = annotation('arrow', ...
            'Color', clr,...
            'headStyle','cback1','HeadLength',headLength/10*norm([dx(ii) dy(ii)]),'HeadWidth',headWidth/10*norm([dx(ii) dy(ii)]));
        %     set(ah,'parent',gca);
        set(ah,'position',[x1(ii) y1(ii) x2(ii)-x1(ii) y2(ii)-y1(ii)]);
        
    end
else
    hkid = get(qh, 'children');
    X = get(hkid(1), 'XData');
    Y = get(hkid(1), 'YData');
    
    
    for ii = 1:3:length(X)-1
        % set the headWidth, function of length of arrow
        dx = diff([X(ii) Y(ii); X(ii+1) Y(ii+1)]);
        ah = annotation('arrow', ...
            'HeadLength', headLength*norm(dx), ...
            'HeadWidth', headWidth*norm(dx), ...
            varargin{:});
        set(ah, 'parent', gca);
        set(ah, 'position', [X(ii) Y(ii) dx]);
    end
end

end
