function Tinv = invertTransform(T, type, n, support)


order = 0;
ct = 0;
while ct < size(T, 1)
    order = order + 1;
    ct = (order + 1) * (order + 2) / 2;
end

if nargin < 2
    type = 'taylor';    
end

if nargin < 3 || isempty(n)
    n = 32;
end

if nargin < 4 || isempty(support)
    support = [-1 1];
end

if ischar(type)
    if strcmpi(type, 'legendre')
        type = @legendreMatrix;
    elseif strcmpi(type, 'taylor')
        type = @taylorMatrix;
    else
        error('Unrecognized basis type');
    end
end

x = linspace(support(1), support(2), n);
y = x;

[X Y] = meshgrid(x, y);

A = type(X(:), Y(:), order);

Xfwd = A*T(:,1);
Yfwd = A*T(:,2);

Afwd = type(Xfwd(:), Yfwd(:), order);

cx = leastSquaresFit(Afwd, X(:));
cy = leastSquaresFit(Afwd, Y(:));

Tinv = cat(2, cx, cy);