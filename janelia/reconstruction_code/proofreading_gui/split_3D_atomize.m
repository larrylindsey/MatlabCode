function split_3D_atomize(body_id, seed_image, seed_depth)
% Segment a designated body based on seeds corresponding to connected
% components in a section
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  11172008  init. code
%

global al cat superpixel_2_seg_map proof data data_splitter_3D

% Create volume in which the 3D splitting is to be performed it has not
% been done yet
if(~isfield('data_splitter_3D', 'volume_3D_split'))
  image = ...
    al{1}(1:data.splitter_3D.scale_factor_sqrt^2:end, ...
    1:data.splitter_3D.scale_factor_sqrt^2:end);
  data_splitter_3D.volume_3D_split = uint8(zeros(size(image,1), size(image,2), length(al)));
  for i = 1:length(al)
    image = al{i}(1:data.splitter_3D.scale_factor_sqrt:end, ...
      1:data.splitter_3D.scale_factor_sqrt:end);
    image_f = 1 - ...
      medfilt2(histeq(double(image)/255), ...
      [data.splitter_3D.med_filt_size data.splitter_3D.med_filt_size]);
    data_splitter_3D.volume_3D_split(:,:,i) = ...
      uint8(255*image_f(1:data.splitter_3D.scale_factor_sqrt:end, ...
      1:data.splitter_3D.scale_factor_sqrt:end));
  end
end


%
% Routine for computing superpixel to new body assignments
%
% get the slice the click was made in
seed_depth = uint32(seed_depth-1); % C indexing

mode = 2; % for additional functionality. Currently redundant

% pmap = uint32(proof.pmap);
label_new = int32(zeros(size(data_splitter_3D.volume_3D_split)));
gui_3D_split_watershed(mode, data_splitter_3D.volume_3D_split, ...
  data.splitter_3D.scale_factor_sqrt^2, cat, superpixel_2_seg_map, proof.pmap, ...
  body_id, seed_image, seed_depth, label_new);

new_labels = unique(label_new(label_new>0));
new_label_mapping = zeros(max(new_labels),1);
new_label_mapping(new_labels) = (1:length(new_labels)) + double(max(proof.pmap));

undo_data = {};
for z = 1:length(al)
  % get membership of superpixels in the new 3D bodies
  superpixel_to_body_map = assign_superpixel_from_3D_split(cat, label_new, ...
    data.splitter_3D.scale_factor_sqrt^2, double(z-1))';
  if(isempty(superpixel_to_body_map))
    continue;
  end
  superpixel_to_body_map = superpixel_to_body_map(min(superpixel_to_body_map,[],2)>0,:);
  superpixel_to_body_map(:,2) = new_label_mapping(superpixel_to_body_map(:,2));
  % find 2D segments whose superpixels belong to multiple bodies. Split
  % these.
  segment_to_body_map = [];
  segment_to_body_map(:,1) = superpixel_2_seg_map{z}(1+superpixel_to_body_map(:,1));
  segment_to_body_map(:,2) = superpixel_to_body_map(:,2);
  segment_to_body_map_unique = unique(segment_to_body_map, 'rows');
  for i = 1:size(segment_to_body_map_unique,1)-1
    if(segment_to_body_map_unique(i,1)==segment_to_body_map_unique(i+1,1))
      % this segment is mapped to multiple bodies so take this subset of
      % superpixels and map them to a new segment which is mapped to the
      % new body.
      sp_set_0 = superpixel_to_body_map(...
        superpixel_to_body_map(:,2)==segment_to_body_map_unique(i,2), 1);
      sp_set = ...
        sp_set_0(superpixel_2_seg_map{z}(1+sp_set_0)==segment_to_body_map_unique(i,1));
      new_segment_id = length(proof.pmap);
      
      undo_data{end+1} = [repmat(z, [length(sp_set),1]), ...
        1+sp_set, reshape(superpixel_2_seg_map{z}(1+sp_set), [length(sp_set),1])];
      superpixel_2_seg_map{z}(1+sp_set) = new_segment_id;
      proof.pmap(1+new_segment_id) = uint32(segment_to_body_map_unique(i,2));
      proof.tmap(1+new_segment_id) = 1; % unknown type
      proof.ttime(1+new_segment_id,:) = 0;
    else
      % the entire (remaining) segment is mapped to a unique body. Just
      % update this.
      undo_data{end+1} = ...
        [1+segment_to_body_map_unique(i,1), proof.pmap(1+segment_to_body_map_unique(i,1))];
      proof.pmap(1+segment_to_body_map_unique(i,1)) = uint32(segment_to_body_map_unique(i,2));
    end
  end
  % the entire (remaining) segment is mapped to a unique body. Just
  % update this.
  undo_data{end+1} = ...
    [1+segment_to_body_map_unique(end,1), proof.pmap(1+segment_to_body_map_unique(end,1))];
  proof.pmap(1+segment_to_body_map_unique(end,1)) = uint32(segment_to_body_map_unique(end,2));
end

data.mmax = double(max(proof.pmap));

data.undo(end+1).action=18;
data.undo(end).slc=seed_depth+1;
data.undo(end).data=undo_data;
data.undo(end).proof=proof.first;

return
end
