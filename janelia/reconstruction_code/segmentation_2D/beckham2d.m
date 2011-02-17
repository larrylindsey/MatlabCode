function pv = beckham2d(im, ts, invert, areat)
    if nargin < 2,
        ts = 100:200;
    end
    if nargin < 3,
        invert = 1;
    end
    if nargin < 4,
        areat = 0;
    end
    
    if invert,
        im = 255-im;
    end
    im = im2double(im);
    ws = zeros(size(im), 'uint32');
    for z=1:size(im,3),
        ws(:,:,z) = watershed(im(:,:,z), 4);
    end
    ws = double(ws);
    pv = zeros(size(im), 'uint8');
    for t=ts,
        l = zeros(size(im), 'uint8');
        for z=1:size(im,3),
            l(:,:,z) = ~remove_merged_boundaries_2D(uint32( ...
                compute_segmentation_from_superpixels_with_min_area_c( ...
                                    ws(:,:,z), im(:,:,z), t/255, areat)));
        end
        pv(logical(l)) = t;
    end
end