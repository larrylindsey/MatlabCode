function [imsplit mask] = splitImageByFlaws(imname, method)

if ischar(imname)
    im = im2single(imreadgray(imname));
else
    im = imname;
end

im_map = 4 * (im - im.*im);

if method == 1
    im_map([1 end],:) = 0;
    im_map(:,[1 end]) = 0;
    im_map_bwl = bwlabel(im_map == 0, 4);
    im_mask = not(im_map_bwl == im_map_bwl(1));
    
    im_chunk_cc = bwconncomp(im_mask);
    goodlabels = false(1, numel(im_chunk_cc.PixelIdxList));
    nth = numel(im) / 16;
    
    for ii = 1:numel(goodlabels)
        goodlabels(ii) = numel(im_chunk_cc.PixelIdxList{ii}) >= nth;
    end
    
    goodlabels = find(goodlabels);
    n = numel(goodlabels);
    
    imsplit = zeros([size(im) n]);
    mask = false([size(im) n]);
    for ii = 1:n
        iml = zeros(size(im));        
        maskl = false(size(im));
        sel = im_chunk_cc.PixelIdxList{goodlabels(ii)};
        iml(sel) = im(sel);
        maskl(sel) = true;
        imsplit(:,:,ii) = iml;
        mask(:,:,ii) = maskl;
    end
elseif method == 2
    vSeamCost = calculateCostmap(im_map) / size(im_map, 1);
    hSeamCost = calculateCostmap(im_map') / size(im_map, 2);
    
    %poly0 = [0 0; 0 size(im, 2); size(im); size(im, 1) 0];
    
    if min(hSeamCost(end,:)) < min(vSeamCost(end,:))        
        seam = calculate_seam(im_map', hSeamCost);
        pseam = cat(1, 1:size(im_map, 2), seam);
        
        poly = cat(2, [size(im,2) 1; 1 1]', pseam);
    else
        seam = calculate_seam(im_map, vSeamCost);
        pseam = cat(1, seam, 1:size(im_map, 1));
               
        poly = cat(2, [1; 1], pseam, [1; size(im,1)]);
    end
    
    seampix = im_map(sub2ind(size(im_map), pseam(2,:), pseam(1,:)));
    
    if median(seampix) * 2 < median(im_map(:))
        mask = poly2mask(poly(1,:), poly(2,:), size(im,1), size(im, 2));
        mask = cat(3, mask, not(mask));
        
        imsplit = cat(3, im, im);
        
        imsplit(mask) = 0;
    else
        mask = true(size(im));
        imsplit = im;
    end
end



