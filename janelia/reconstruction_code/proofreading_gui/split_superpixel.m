function split_superpixel()
% function to split superpixel based on user-input seeds using Graph Cut
% algorithm. Called by gui.m v2.8.3+
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  ~05262008   init code.
%

global cat al data superpixel_2_seg_map proof


min_x = min(data.superpixel_splitter.member_pixel_pos(:,1));
max_x = max(data.superpixel_splitter.member_pixel_pos(:,1));
min_y = min(data.superpixel_splitter.member_pixel_pos(:,2));
max_y = max(data.superpixel_splitter.member_pixel_pos(:,2));

% get the image cropping for the superpixel
image_superpixel = al{data.superpixel_splitter.ik}(min_y:max_y, min_x:max_x);

data.superpixel_splitter.member_pixel_pos(:,1) = data.superpixel_splitter.member_pixel_pos(:,1) - min_x+1;
data.superpixel_splitter.member_pixel_pos(:,2) = data.superpixel_splitter.member_pixel_pos(:,2) - min_y+1;

mask_superpixel = zeros(size(image_superpixel));
mask_superpixel((data.superpixel_splitter.member_pixel_pos(:,1)-1)*size(image_superpixel,1)+data.superpixel_splitter.member_pixel_pos(:,2)) = 1;

% boundary map that would determine the pairwise energy terms
boundary_map = int32(255*medfilt2(histeq(im2double(image_superpixel)), [10 10]).*mask_superpixel);

% user input seeds for label 0 and 1
if(isempty(data.superpixel_splitter.mouse_points_label_0) || isempty(data.superpixel_splitter.mouse_points_label_1))
  return;
end;

data.superpixel_splitter.mouse_points_label_0 = round(data.superpixel_splitter.mouse_points_label_0);
data.superpixel_splitter.mouse_points_label_0(:,1) = data.superpixel_splitter.mouse_points_label_0(:,1) - min_x + 1;
data.superpixel_splitter.mouse_points_label_0(:,2) = data.superpixel_splitter.mouse_points_label_0(:,2) - min_y + 1;
is_valid_x = data.superpixel_splitter.mouse_points_label_0(:,1)>=1 & data.superpixel_splitter.mouse_points_label_0(:,1)<=size(image_superpixel,2);
is_valid_y = data.superpixel_splitter.mouse_points_label_0(:,2)>=1 & data.superpixel_splitter.mouse_points_label_0(:,2)<=size(image_superpixel,1);
data.superpixel_splitter.mouse_points_label_0 = data.superpixel_splitter.mouse_points_label_0(is_valid_x & is_valid_y,:);

data.superpixel_splitter.mouse_points_label_1 = round(data.superpixel_splitter.mouse_points_label_1);
data.superpixel_splitter.mouse_points_label_1(:,1) = data.superpixel_splitter.mouse_points_label_1(:,1) - min_x + 1;
data.superpixel_splitter.mouse_points_label_1(:,2) = data.superpixel_splitter.mouse_points_label_1(:,2) - min_y + 1;
is_valid_x = data.superpixel_splitter.mouse_points_label_1(:,1)>=1 & data.superpixel_splitter.mouse_points_label_1(:,1)<=size(image_superpixel,2);
is_valid_y = data.superpixel_splitter.mouse_points_label_1(:,2)>=1 & data.superpixel_splitter.mouse_points_label_1(:,2)<=size(image_superpixel,1);
data.superpixel_splitter.mouse_points_label_1 = data.superpixel_splitter.mouse_points_label_1(is_valid_x & is_valid_y,:);


% construct the label energy map from the user-input seeds
label_energy_prior_1 = zeros(size(image_superpixel));
label_energy_prior_1((data.superpixel_splitter.mouse_points_label_0(:,1)-1)*size(image_superpixel,1) + ...
  data.superpixel_splitter.mouse_points_label_0(:,2)) = 100;
label_energy_prior_1 = int32(imdilate(label_energy_prior_1, strel('square', 2*data.superpixel_splitter.brush_size)));

label_energy_prior_0 = zeros(size(image_superpixel));
label_energy_prior_0((data.superpixel_splitter.mouse_points_label_1(:,1)-1)*size(image_superpixel,1) + ...
  data.superpixel_splitter.mouse_points_label_1(:,2)) = 100;
label_energy_prior_0 = int32(imdilate(label_energy_prior_0, strel('square', 2*data.superpixel_splitter.brush_size)));

label_energy_prior = int32([reshape(label_energy_prior_0,[1 numel(image_superpixel)]); ...
  reshape(label_energy_prior_1,[1 numel(image_superpixel)])]);


% for debuggin purposes
% figure(1001); imshow(image_superpixel);
% figure(1002); imshow(boundary_map, []);
% figure(1003); imshow(label_energy_prior_1, []);
% figure(1004); imshow(label_energy_prior_0, []);

n_iteration = data.superpixel_splitter.n_iteration;

% perform Graph Cut
label_map = segmentation_2D_graphcut_2(boundary_map, label_energy_prior, n_iteration);
% display_image = repmat(image_superpixel, [1 1 3]); display_image(:,:,3) = 255*label_map;
% figure(1005); imshow(display_image);

% % get all the connected components in the two labelling and assign them to
% % different bodies. The largest component in label_map==0 retains its
% % body-id.
% % First label 0
% label_0_mask = label_map==0;


% modification ids
new_body_id = max(proof.pmap)+1;
new_segment_id = length(proof.pmap);
new_superpixel_id = max(cat{data.superpixel_splitter.ik}(:))+1;
[py, px] = find(mask_superpixel);
py = py + min_y - 1;
px = px + min_x - 1;
label_map_1 = label_map;
label_map(label_map_1==0) = cat{data.superpixel_splitter.ik}(py(1),px(1));
label_map(label_map_1==1) = new_superpixel_id;
area_new_superpixel = nnz(label_map_1);
label_boundary = draw_segment_boundaries_c(label_map);
label_map(label_boundary==1) = 0;
pixel_ids = (px-1)*size(cat{data.superpixel_splitter.ik},1) + py;

% undo buffer
i=length(data.undo);
data.undo(i+1).action=14;
data.undo(i+1).slc=data.superpixel_splitter.ik;
data.undo(i+1).data={pixel_ids, cat{data.superpixel_splitter.ik}(pixel_ids)};
data.undo(i+1).proof=proof.first;

% linkage graph buffers
curr_superpixel_id = cat{data.superpixel_splitter.ik}(pixel_ids(1));
bodies_modified = [proof.pmap(1+superpixel_2_seg_map{data.superpixel_splitter.ik}(1+curr_superpixel_id)) ...
  new_body_id];
for b_m = bodies_modified
  data.links_3D_update.buffer(end+1).body_id = b_m;
  data.links_3D_update.buffer(end).z = data.superpixel_splitter.ik;
end

% update the superpixel map and the superpixel_2_seg_map
cat{data.superpixel_splitter.ik}(pixel_ids) = label_map(mask_superpixel==1);
superpixel_2_seg_map{data.superpixel_splitter.ik}(new_superpixel_id+1) = new_segment_id;
proof.pmap(new_segment_id+1) = new_body_id;
proof.tmap(new_segment_id+1)=1;
proof.ttime(new_segment_id+1,:) = 0;

update_for_new_segment(new_segment_id, new_body_id, area_new_superpixel, data.superpixel_splitter.ik);

% remember in .tmp{20}
if(isempty(data.tmp{20}))
  data.tmp{20}=cell(1,length(cat));
end
if(isempty(data.tmp{20}{data.superpixel_splitter.ik}))
  data.tmp{20}{data.superpixel_splitter.ik}{1}=[];
  data.tmp{20}{data.superpixel_splitter.ik}{2}=[];
end

data.tmp{20}{data.superpixel_splitter.ik}{1} = ...
  [data.tmp{20}{data.superpixel_splitter.ik}{1}(:); pixel_ids];
data.tmp{20}{data.superpixel_splitter.ik}{2} = ...
  [data.tmp{20}{data.superpixel_splitter.ik}{2}(:); double(label_map(mask_superpixel==1))];

return
end
