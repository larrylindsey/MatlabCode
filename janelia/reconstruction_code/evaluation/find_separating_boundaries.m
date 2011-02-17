function boundary_map = find_separating_boundaries(segmentation, body_pairs_list)
% BOUNDARY_MAP = FIND_SEPARATING_BOUNDARIES(SEGMENTATION, BODY_PAIRS_LIST)
% For an n-by-2 list of body labels, find boundaries separating segments 
% BODY_PAIRS_LIST(:,1) and BODY_PAIRS_LIST(:,2) in SEGMENTATION
    boundary_map = zeros(size(segmentation), 'uint8');
    se = zeros([3 3 3]);
    se(:,2,2) = 1;
    se(2,:,2) = 1;
    se(2,2,:) = 1;
    for i=1:size(body_pairs_list, 1),
        l1 = body_pairs_list(i, 1);
        l2 = body_pairs_list(i, 2);
        boundary = logical(imdilate(segmentation == l1, se) ...
                        .* imdilate(segmentation == l2, se) ...
                        .* (segmentation == 0) ...
                        );
        boundary_map(boundary) = 1;
    end
end