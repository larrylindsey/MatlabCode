function label_map = watershed_seeded(image, seed_map)
% label_map = watershed_seeded(image, seed_map)
% Perform seeded 2D watershed using imimposemin. Segment labels are
% assigned according to seed_map.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  04202009  init. code
%

seed_mask = seed_map>0;
image_seeded = imimposemin(image, seed_mask);
seg_ws = watershed(image_seeded);

label_pair = [seg_ws(:), seed_map(:)];
label_pair_u = unique(label_pair, 'rows');
label_pair_u = label_pair_u(min(label_pair_u, [], 2)>0, :);

seg_ws_mapping = zeros(max(seg_ws(:))+1, 1);
seg_ws_mapping(label_pair_u(:,1)+1) = label_pair_u(:,2);
label_map = seg_ws_mapping(seg_ws+1);
return
end
