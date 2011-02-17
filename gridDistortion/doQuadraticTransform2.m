function U = doQuadraticTransform2(X, T, select)

if nargin < 3
    select = false;
end

%X in P by 2
%T in 6 by 2 or 12 by 2

if isa(X, 'logical')
    if X
        select = 'true';
    else
        select = 'false';
    end
    U = eval(['@(X, T) doQuadraticTransform2(X, T, ' select ');']);
else    
    %addIndex = select * 6;

    if isstruct(T)
        T = T.tdata;
    end

    if select
        T = T(7:12, :);
    else
        T = T(1:6, :);
    end

    
    x2 = X(:, 1) .^ 2;
    y2 = X(:, 2) .^ 2;
    xy = X(:, 1) .* X(:, 2);
    x1 = X(:, 1);
    y1 = X(:, 2);
    o = ones(size(X, 1), 1);
    
    A = cat(2, x2, xy, y2, x1, y1, o);
   
    %A in P by 6, T in 6 by 2
    
    U = A * T;
    
    %U in P by 2
end