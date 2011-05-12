function [contours iout] = maskStackToContours(stack, x, y, iin)


sz = size(stack);
n = sz(3);
sz = sz(1:2);

if nargin < 4
   iin = 1:n;
end


contours = {};

x0 = x(1);
xr = x(2) - x(1);
y0 = y(1);
yr = y(2) - y(1);

iout = [];

for ii = 1:n
    currmask = logical(stack(:,:,ii));
    currmask = imclose(currmask, strel('disk', 1));
    
    P = bwboundaries(currmask, 8);
    
    for ip = 1:numel(P)
        P{ip} = P{ip} ./ repmat(sz, [size(P{ip}, 1) 1]);
        P{ip}(:,2) = P{ip}(:,2) * xr + x0;
        P{ip}(:,1) = (1 - P{ip}(:,1)) * yr + y0;
        P{ip} = P{ip}(:,[2 1]);
    end
    
    contours = cat(1, contours, P);
    
    iout = cat(1, iout, repmat(iin(ii), [numel(P) 1]));
end

% for ii = 1:numel(contours)
%     p = contours{ii};
%     p = p ./ repmat(sz, [size(p, 1) 1]);
%     p(:,2) = p(:,2) * xr + x0;
%     p(:,1) = p(:,1) * yr + y0;
%     p = p(:,[2 1]);
%     contours{ii} = p;
% end

end