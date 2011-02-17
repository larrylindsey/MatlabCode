function Tout = xformFromReconstruct(Tin, mag)

if size(Tin, 2) > size(Tin, 1)
    Tin = Tin.';
end

Tin = Tin([1 2 3 5 4 6], :);
Tout = cat(1, Tin, zeros(4, 2));

Tout = flipud(fliplr(Tout));

if nargin > 1
    magmat = repmat(mag, size(Tout));
    powmat = repmat([3 3 3 3 2 2 2 1 1 0]', [1 2]);
    magmat = magmat .^ (powmat - 1);

    Tout = Tout .* magmat;
end