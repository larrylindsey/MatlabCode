function Ttrans = transformTransform(T, transfuncXY, type, n, support)


if nargin < 5
    support = [-1 1];
    if nargin < 4
        n = 32;
        if nargin < 3
            type = 'taylor';
        end
    end
end


if ischar(type)
    if strcmpi(type, 'legendre')
        type = @legendreMatrix;
    elseif strcmpi(type, 'taylor')
        type = @taylorMatrix;
    else
        error('Unrecognized basis type');
    end
elseif ~isa(type, 'function_handle')
    error('Unrecognized basis type');
end


if isa(transfuncXY, 'function_handle')
    transfuncXY = {transfuncXY, transfuncXY};
elseif isa(transfuncXY, 'cell')
    if ~isa(transfuncXY{1}, 'function_handle')
        error('Unrecognized transform type');
    end
else
    error('Unrecognized transofmr type');
end

x = linspace(support(1), support(2), n);
y = x;

[X Y] = meshgrid(x, y);

A = type(X(:), Y(:));

Xfwd = A*T(:,1);
Yfwd = A*T(:,2);

Xfwd = transfuncXY{1}(Xfwd);
Yfwd = transfuncXY{2}(Yfwd);

cx = leastSquaresFit(A, Xfwd(:));
cy = leastSquaresFit(A, Yfwd(:));

Ttrans = cat(2, cx, cy);