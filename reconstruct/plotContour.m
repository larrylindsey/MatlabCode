function handles = plotContour(varargin)
%name, secdoc, fig

switch numel(varargin)
    case 1
        contour = varargin{1};
        fig = figure;
        doview = true;
    case 2
        if isstruct(varargin{2})
            contour = extractContour(varargin{2}, varargin{1});
            fig = figure;
            doview = true;
        elseif isnumeric(varargin{2})
            contour = varargin{1};
            doview = false;
            fig = varargin{2};
        else
            usageError;
        end
    case 3
        contour = extractContour(varargin{2}, varargin{1});
        doview = false;
        fig = varargin{3};
    otherwise
        usageError;       
end

%
% if nargin < 3
%     fig = figure;    
%     doview = true;        
% else
%     doview = false;
% end
% 
% 
% contour = extractContour(secdoc, name);
%dz = secdoc(1).section.thickness;

figure(fig); hold on;

hh = zeros(1, numel(contour));
for i_c = 1:numel(contour)
    x = contour(i_c).transPoints(:,1);
    y = contour(i_c).transPoints(:,2);
    if contour(i_c).closed
        x(end + 1) = x(1);
        y(end + 1) = y(1);
    end
    z = repmat(contour(i_c).z, size(x));
    

    hh(i_c) = plot3(x, y, z, 'Color', contour(i_c).border);
end

axis equal;
grid on;

if doview
    view(45, 45);
end

if nargout > 0
    handles = hh;
end

end

function usageError
error(['Usage: plotContour(name, secdoc [,fig]),'...
                ' or\n       plotContour(contour [,fig])']);
end