function tr = transposeTransform(tr)

o = tr.order;

n = (o + 1) * (o + 2) / 2;

jj = 1;
th = 1;
k = 1;
next_th = 2;

sel = zeros(1, n);

for ii = 1:n
    sel(ii) = k;
    k = k - 1;
    if k < th
        jj = jj + 1;
        k = jj * (jj + 1) / 2;
        th = next_th;
        next_th = k + 1;
    end
end

tr.T = tr.T(sel,[2 1]);
tr.Tinv = tr.Tinv(sel,[2 1]);

if isfield(tr, 'data') && isfield(tr.data, 'u') && isfield(tr.data, 'v')
    tmp = tr.data.u;
    tr.data.u = tr.data.v;
    tr.data.v = tmp;
end
