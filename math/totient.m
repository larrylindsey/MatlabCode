function t = totient(x)

t = zeros(size(x));

for it = 1:numel(x)
    g = gcd(1:(x(it)-1), x(it));
    t(it) = sum(g == 1);
end
