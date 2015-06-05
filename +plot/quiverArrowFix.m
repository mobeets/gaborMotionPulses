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
