function imGD = quickImDistort(im, K, fill)

sz = size(im);
sz = sz(1:2);

[X Y] = meshgrid(1:sz(1), 1:sz(2));

[xp yp] = distort(X(:), Y(:), K, sz(1) / 2, sz(2) / 2);
xpi = round(xp);
ypi = round(yp);
baduns = logical(xpi < 1);
baduns = or(baduns, logical(ypi < 1));
baduns = or(baduns, logical(xpi > sz(1)));
baduns = or(baduns, logical(ypi > sz(2)));
xpi(baduns) = 1;
ypi(baduns) = 1;
imGD = zeros(size(xpi));
for ii = 1:numel(xpi)
imGD(ii) = im(xpi(ii), ypi(ii));
end
imGD(baduns) = fill;
imGD = reshape(imGD, sz);
