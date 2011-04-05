function [xu yu] = functionify(x, y)
% [x y] = functionify(x, y)
% Removes non-unique assignment for values x
% y's returned for non-unique x's are taken to be the mean.

it = 1:numel(x); %all indices

[x isort] = sort(x);
y = y(isort);

[xu iu] = unique(x);
yu = y(iu);
inu = setdiff(it, iu);

xnu = x(inu);

for ix = 1:numel(xnu)
    insel = x == xnu(ix);
    outsel = xu == xnu(ix);
    yu(outsel) = mean(y(insel));
end