function Aout = taylorMat(x, y, order)

if nargin < 3
    order = 3;
end

x = x(:);
y = y(:);

o = (order + 1) * (order + 2) / 2;

Aout = zeros(numel(x), o);

ii = 1;
for i_order = 0:order
    for n = 0:i_order
        x_order = i_order - n;
        y_order = n;
        Aout(:,ii) = (x.^x_order) .* (y.^y_order);
        ii = ii + 1;
    end    
end

% 
% x1 = x(:);
% y1 = y(:);
% 
% %x1, y1 in [n 1]
% 
% x2 = x1.*x1;
% xy = x1.*y1;
% y2 = y1.*y1;
% 
% x3 = x2.*x1;
% x2y = x2.*y1;
% xy2 = x1.*y2;
% y3 = y2.*y1;
% 
% o = ones(size(x1));
% 
% Aout = cat(2, o, x1, y1, x2, xy, y2, x3, x2y, xy2, y3);