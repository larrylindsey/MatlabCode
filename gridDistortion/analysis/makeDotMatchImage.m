function imdot = makeDotMatchImage(r, im, rc_found, rc_grid, grid_model)

imdotR = im2double(im);
if size(imdotR, 3) > 1
    imdotR = rgb2gray(imdotR);
end
imdotG = imdotR;
imdotB = imdotR;

n = size(rc_found, 1);

imDotFound = full(sparse(round(rc_found(:,1)), round(rc_found(:,2)),...
    ones(n, 1), size(im,1), size(im,2)));
imDotFound = logical(imDotFound);
imDotFound = imdilate(imDotFound, strel('disk', r));

if nargin > 3

    rc_match = linearAlign(rc_grid * grid_model, rc_found);
    imDotMatch = full(sparse(round(rc_match(:,1)), round(rc_match(:,2)),...\
        ones(n, 1), size(im,1), size(im,2)));
    imDotMatch = logical(imDotMatch);
    imDotMatch = imdilate(imDotMatch, strel('disk', r));
else
    imDotMatch = false(size(imDotFound));
end


imdotR(imDotMatch) = 0;
imdotG(or(imDotMatch, imDotFound)) = 0;
imdotB(imDotFound) = 0;

imdotR(imDotFound) = 1;
imdotB(imDotMatch) = 1;

imdot = cat(3, imdotR, imdotG, imdotB);