function tr = identityTransform(type, order, d, n, lim)


if nargin < 5
    lim = [-1; 1];
    if nargin < 4
        n = 1024;
        if nargin < 3
            d = 2;
            if nargin < 2
                order = 1;
                if nargin < 1
                    type = @taylorMat;
                end
            end
        end
    end
end

if size(lim, 2) == 1
    lim = repmat(lim, [1 d]);
end
if numel(n) == 1
    n = repmat(round(n^(1/d)), [d 1]);
end

t = cell(1,d);

for ii = 1:d
    t{ii} = linspace(lim(1,ii), lim(2,ii), n(ii));
end

rc = gridRC(t{:});

tr = fitTransform(rc, rc, order, type);

