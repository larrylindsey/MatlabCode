function [im_out sz_out] = greedyResize(im, reduction, efunc)

if nargin < 3
    efunc = @defaultEnergy;
end

sz = size(im);
sz = sz(1:2);

newsz = sz - reduction;

sz_out = sz;

if isinteger(im)
    im = double(im) / 255;
end

im_out = im;

while prod(sz(1:2) - newsz) > 0
    energy = efunc(im_out);
    costmap_c = calculateCostmap(energy);
    costmap_r = calculateCostmap(energy');
    minc = min(costmap_c(end,:));
    minr = min(costmap_r(end,:));
    if minc < minr
        im_out = reduceWidth(im_out, 1, efunc, costmap_c);
    else
        im_out = reduceHeight(im_out, 1, efunc, costmap_r);
    end
    sz = size(im_out);
    sz_out = cat(1, sz_out, sz(1:2));
end

if sz(1) > newsz(1)
    [im_out szhist] = reduceHeight(im_out, sz(1) - newsz(1), efunc);
    sz_out = cat(1, sz_out, szhist);
elseif sz(2) > newsz(2)
    [im_out szhist] = reduceWidth(im_out, sz(2) - newsz(2), efunc);    
    sz_out = cat(1, sz_out, szhist);
end



end