function [cost, p] = getMinImagePath(im, rc1, rc2)

% Intended for use with small images

% A = sparse(numel(im), numel(im));
ii = reshape(1:numel(im), size(im));
iisz = size(ii);
dv = [-1 -1; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1];
dn = sqrt(sum(dv.^2, 2));

is = zeros([size(im, 1), size(im,2), size(dv,1)]);
js = is;
as = is;

for r = 1:size(im,1)
    for c = 1:size(im,2)
        for d = 1:size(dv,1)
            ro = r + dv(d,1);
            co = c + dv(d,2);
            if rcok(ro, co, iisz)
                is(r, c, d) = ii(r,c);
                js(r, c, d) = ii(ro,co);
                as(r, c, d) = dn(d) * (im(r,c) + im(ro,co))/2;
            end
        end
    end
end

sel = is(:) > 0;
is = is(sel);
js = js(sel);
as = as(sel);

A = sparse(is, js, as, numel(im), numel(im));

[sp_cost, sp_pred] = thirdpartyfunction('shortest_paths',...
    A, ii(rc1(1), rc1(2)));
cost = sp_cost(ii(rc2(1), rc2(2)));

if nargout > 1    
    p0 = ii(rc2(1), rc2(2));
    p = [];
    
    while p0 > 0
        p = [p p0]; %#ok<AGROW>
        p0 = sp_pred(p0);
    end
end
end

function ok = rcok(r, c, sz)
rok = and(r > 0, r <= sz(1));
cok = and(c > 0, c <= sz(2));
ok = and(rok, cok);
end
