function [label, prob] = makeRayLabel(probmap, seedval)

sz = size(probmap);

wslabel = zeros(sz, 'uint16');
lmax = zeros(sz(3), 1);

parfor ii = 1:sz(3)
    pmap = probmap(:,:,ii);    
%     pmap = imimposemin(pmap, pmap <= seedval);
    pmap = medfilt2(pmap, [3 3]);
    lii = watershed(pmap, 4);
    wslabel(:,:,ii) = lii;
    lmax(ii) = max(lii(:));
end

csum = cumsum(lmax);

parfor ii = 2:sz(3)
    lii = wslabel(:,:,ii);
    zmask = lii == 0;
    lii = lii + csum(ii - 1);
    lii(zmask) = 0;
    wslabel(:,:,ii) = lii;
end

label = zeros([sz(1:2), 2 * sz(3) - 1], class(wslabel));
label(:,:,1:2:end) = wslabel;

prob = zeros(size(label), class(probmap));

prob(:,:,1:2:end) = probmap;

end
