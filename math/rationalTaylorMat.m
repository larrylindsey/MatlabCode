function Aout = rationalTaylorMat(xi, yi, to, order)

if numel(order) == 1
    sprintf(['Singleton order.', ...
        ' Assuming equal nominator/denominator order = %d\n'], order);
    order = [order order];
end

P = taylorMat(xi, yi, order(1));
Q = -taylorMat(xi, yi, order(2));

Q(:,1) = [];

%Aout = zeros(size(P) + size(Q));

to = reshape(to, [numel(to) 1]);

Q = Q .* repmat(to, [1 size(Q, 2)]);


%Aout(1:size(P,1), 1:size(P,2)) = P;
%Aout(size(P,1) + (1:size(Q,1)), size(P,2) + (1:size(Q,2))) = Q;
Aout = cat(2, P, Q);
