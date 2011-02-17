function m = makeLineModel(inPoints)
A = inPoints;
c = ones(size(A, 1), 1);

At = A';
AtA = At * A;
[U S V] = svd(AtA);
AtAInv = U * inv(S) * V';
m = AtAInv * At * c;