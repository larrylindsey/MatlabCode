function pn = generateLegendre(n, x, pn_1, pn_2)

pn = ((2 * n - 1) * x .* pn_1 - (n - 1) * pn_2) / n;