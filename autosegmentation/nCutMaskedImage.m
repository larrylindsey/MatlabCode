function [W im mask offset] = nCutMaskedImage(mask, im, nbd, sig)
global sigI;
global sigX;

% Prepare inputs
if nargin < 2
    im = im2double(mask);
    nbd = 3;
    sig = [1 1];
end

if size(im,3) > 1
    im = rgb2gray(im);
end

if numel(nbd) == 1
    [x y] = meshgrid(-nbd:nbd, -nbd:nbd);
    sel = (x.^2 + y.^2) < nbd^2;
    sel = and(sel, or(x ~= 0, y~=0));    
    x = x(sel);
    y = y(sel);
    nbd = cat(2, x(:), y(:));
end

if islogical(nbd)
    [x y] = find(nbd);
    nbd = cat(2, x(:), y(:));
    nbd = nbd - repmat(round(mean(nbd, 1)), [size(nbd, 1) 1]);
    sel = prod(nbd, 2) ~= 0;
    nbd = nbd(sel,:);
end

sigI = sig(1);
sigX = sig(2);


szMask = size(mask);
szIm = size(im);

if numel(szMask) ~= numel(szIm) || any(szMask ~= szIm)
    error('Image and mask must be of the same size');
end

im = im2double(im);

% Crop out unneccesary pixels.
[mask im offset] = doCrop(mask, im);

% Calculate indices
iNodes = calculateIndices(mask);

%Calculate weight-connection matrix
W = calculateW(mask, im, iNodes, nbd);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask im offset] = doCrop(mask, im)
rInd = sum(mask, 2) > 0;
cInd = sum(mask, 1) > 0;

mask = mask(rInd, cInd);
im = im(rInd, cInd);

offset(1) = find(rInd, 1, 'first');
offset(2) = find(cInd, 1, 'first');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iNodes = calculateIndices(mask)
rcs = cumsum(mask, 1);
ccs = circshift(cumsum(rcs(end,:)), [1 1]);
ccs(1) = 0;

ccs = repmat(ccs, [size(rcs,1), 1]);
iNodes = (rcs + ccs) .* mask;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = calculateW(mask, im, iNodes, nbd)
%nbd in [x 2]

imsz = size(im);
imsel = im(mask);
npix = sum(mask(:));
nnbd = size(nbd, 1);
nbdsz = max(abs(nbd), [], 1);

% Index Neighbors
%inbd = zeros([imsz nnbd]);



inbd = zeros([imsz(1) + nbdsz(1) imsz(2) + nbdsz(2) nnbd]);
inbd(1:imsz(1), 1:imsz(2), :) = repmat(iNodes, [1 1 nnbd]);

for is = 1:nnbd
    inbd(:,:,is) = circshift(inbd(:,:,is), nbd(is,:));
end
inbd = inbd(1:imsz(1), 1:imsz(2), :);

spi = cell(1,nnbd);
spj = cell(1,nnbd);
spv = cell(1,nnbd);

for il = 1:nnbd
    nbdslice = inbd(:,:,il);
    
    sel = false(numel(im), 1);
    sel(:) = and(nbdslice > 0, mask);
    
    spv{il} = weightFun(imsel(iNodes(sel)), imsel(nbdslice(sel))) * ...
        distFun(nbd(il,:));
    
    spi{il} = iNodes(sel);
    spj{il} = nbdslice(sel);
    
end

spi = cat(1, spi{:});
spj = cat(1, spj{:});
spv = cat(1, spv{:});

W = sparse(spi, spj, spv, npix, npix);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = weightFun(a, b)
global sigI;

%w = exp(- (a - b).^2 / sigI);
w = a .* b / sigI;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = distFun(v)
global sigX;

w = exp(-sum(v.^2) / sigX);
end