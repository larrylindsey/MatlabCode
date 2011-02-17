function S = cubicbsplineinterp(path, n)

q = size(path, 1) - 4;
t = linspace(0, 1, n + 1);
t = t(1:n);
T = [t.^3; t.^2; t; ones(size(t))]';
S = [];

M = [-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0] / 6;

for jj = 1:q
    p = path((0:3) + jj, :);
    S = cat(1, S, T * M * p);
end