function pred_vol = tsegments2beckham(ws_filename, ...
                                        amb_file_pattern, threshold_list, d)
% PRED_VOL = TSEGMENTS2BECKHAM(WS_FILENAME, AMB_FILE_PATTERN 
%                                                      [, THRESHOLD_LIST, D])
% Generate a volume of Beckham values (Boundary Elimination Criterion At
% Merge) from a group of thresholded segmentations. The function assumes
% that the segmentations have 0-labeled voxels separating different labels.
% The length of the (optional) threshold_list should equal the number of
% files matching amb_file_pattern (and threshold(i) should correspond to the
% ith file as given by dir(). The optional d argument specifies whether the
% data are 3D or a stack of (unlinked) 2D sections.

    if nargin == 2,
        threshold_list = 0;
    end
    if nargin < 4,
        d = 3;
    end
    
    ws = fread_raw_array_mex(expandtilde(ws_filename));
    
    dir_delimiters = find(amb_file_pattern == '/');
    if numel(dir_delimiters) > 0,
        directory = amb_file_pattern(1:dir_delimiters(end));
    else
        directory = './';
    end
    a_files = dir(amb_file_pattern);
    
    if ~threshold_list,
        threshold_list = 1:numel(a_files);
    end
    
    h = waitbar(0, 'Building Beckham map...');
    
    pred_vol = uint8(~read_segmentation_from_raw(ws, ...
                                            [directory a_files(1).name], d));
    
    waitbar(1/numel(a_files), h);
    
    for i=2:numel(a_files),
        cur_pred_vol = uint8(~read_segmentation_from_raw(ws, ...
                                            [directory a_files(i).name], d));
        cur_thresh = threshold_list(i);
        update_ids = logical((cur_pred_vol ~= 0).*(pred_vol < cur_thresh));
        pred_vol(update_ids) = cur_thresh;
        waitbar(i/numel(a_files), h);
    end
    delete(h);
end