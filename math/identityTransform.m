function tr = identityTransform(order, type, n, u, v)


if nargin < 5
    v = [-1 1];
    if nargin < 4
        u = [-1 1];
        if nargin < 3
            n = 64;            
            if nargin < 2
                type = @legendreMat;
            end
        end
    end
end

rc = gridRC(linspace(u(1), u(2), n), linspace(u(1), u(2), n));

data.n = n;
data.u = u;
data.v = v;

tr = regressionTransform(rc, rc, order, type, data);
tr.T(tr.T < 1) = 0;
tr.Tinv(tr.Tinv < 1) = 0;

