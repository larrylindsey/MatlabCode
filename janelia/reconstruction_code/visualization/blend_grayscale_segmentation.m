function imout = blend_grayscale_segmentation(imin, seg, cmap, alpha)
% IMOUT = BLEND_GRAYSCALE_SEGMENTATION(IMIN, SEG [, CMAP, ALPHA])
% Combine grayscale image with a segmentation map into a single color
% image stack. 
% INPUT:
%   - IMIN: the input grayscale image stack (size: [m n p])
%   - SEG: the input segmentation (size: [m n p])
%   - CMAP: (optional) a colormap to use instead of jet.
%   - ALPHA: (optional) the opacity of the segmentation map. Default is 1.
% OUTPUT:
%   - IMOUT: the output RGB color stack (size: [m n 3 p])
    if nargin < 3,
        cmap = [];
    end
    if nargin < 4,
        alpha = 1;
    end
    if strcmp(class(imin), 'uint8'),
        imin = double(imin)/255;
    end
    imin = reshape(imin, [size(imin,1), size(imin,2), 1, size(imin,3)]);
    seg = reshape(seg, [size(seg,1), size(seg,2), 1, size(seg,3)]);
    red = imin;
    green = imin;
    blue = imin;
    voxels_to_update = logical(seg);
    levels = sort(unique(seg(voxels_to_update)));
    if isempty(cmap),
        cmap = jet(numel(levels));
        if size(cmap,1) == 1,
            % If we are only using a single threshold, make the corresponding
            % pixels red
            cmap = [1 0 0];
        end
    end
    h = waitbar(0, 'Computing blended stack...');
    for i=1:numel(levels),
        voxels_to_update = seg == levels(i);
        red(voxels_to_update) = (1-alpha)*imin(voxels_to_update) + ...
                                                            alpha*cmap(i,1);
        green(voxels_to_update) = (1-alpha)*imin(voxels_to_update) + ...
                                                            alpha*cmap(i,2);
        blue(voxels_to_update) = (1-alpha)*imin(voxels_to_update) + ...
                                                            alpha*cmap(i,3);
        waitbar(i/numel(levels), h);
    end
    imout = cat(3, red, green, blue);
    delete(h);
end