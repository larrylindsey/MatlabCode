function segmented_vol = read_segmentation_from_raw(ws, amb, d)
% SEGMENTED_VOL = READ_SEGMENTATION_FROM_RAW(WS, AMB [, D]) 
% Read a segmented volume from the pixel-to-superpixel watershed map and
% the superpixel-to-body (usually AMB) map. These can be either the arrays
% themselves or the filenames. The optional argument "d" can be either 2 or
% 3 and is used only to remove merged boundaries either in 2D or 3D.

    if nargin == 2,
        d = 3;
    end
    
    if ischar(ws),
        ws_fn = ws;
        ws = fread_raw_array_mex(expandtilde(ws));
        if size(ws) == 0,
            error('ws file reading failed (size=0), filename:\n%s', ws_fn);
        end
    end
    if ischar(amb),
        amb_fn = amb;
        amb = fread_raw_array_mex(expandtilde(amb));
        if size(amb) == 0,
            error('amb file reading failed (size=0), filename:\n%s', amb_fn);
        end        
    end
    
    m = zeros(max(amb(1,:))+1,1);
    m(1+amb(1,:)) = amb(2,:);
    
    segmented_vol_t = uint32(m(1+ws));
    if d == 3,
        segmented_vol_t = remove_merged_boundaries_3D(segmented_vol_t);
    elseif d == 2,
        for i=1:size(segmented_vol_t, 3),
            segmented_vol_t(:,:,i) = ...
                            remove_merged_boundaries_2D(segmented_vol_t(:,:,i));
        end
    else
        error('Number of dimensions d should be 2 or 3 (%d given).\n', d);
    end
    
    % The array needs to be transposed because Matlab uses column-major
    % storage instead of the C standard row-major.
    segmented_vol = zeros([size(segmented_vol_t, 2), ...
        size(segmented_vol_t, 1), size(segmented_vol_t, 3)], ...
        class(segmented_vol_t));
    for i=1:size(segmented_vol, 3),
        segmented_vol(:,:,i) = segmented_vol_t(:,:,i)';
    end
end