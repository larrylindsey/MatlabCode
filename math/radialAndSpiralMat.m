function [A1, A2] = radialAndSpiralMat(x, order)

if size(x, 2) ~= 2
    error('radialAndSpiralMat supports only 2 dimensions, rather than %d', ...
        size(x, 2));
end

r = sqrt(sum(x.^2, 2));

P = taylorMatND(r, order);
S = taylorMatND(r, 2);

l = size(P,2);

P1 = repmat(x(:,1), [1 l]) .* P;
P2 = repmat(x(:,2), [1 l]) .* P;
S1 = repmat(x(:,2), [1 3]) .* S;
S2 = repmat(x(:,1), [1 3]) .* S;


A1 = cat(2, P1, S1);
A2 = cat(2, P2, S2);

end
