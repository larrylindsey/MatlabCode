function generate_cat_superpixel_2_seg_map_multi_tile(config)
% generate_cat_multi_tile(config)
% generate cat - multiple tiles per section (trakEM xml)
%
% cat{z}'s are the superpixel maps being used in the reconstruction
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified for multiple tiles per section (trakEM xml).
%                 Including generation of superpixel_2_seg_map.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Constructing cat and superpixel_2_seg_map ..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stack_config = config.stack;
seg_config = config.superpixel_2_seg(end);
superpixel_config = config.superpixel(end);
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

sift_dir = [get_region_dir(config), config.SIFT.dir];
align_seg_dir = get_align_seg_dir(config);

align_roi_file_name = ['./', stack_config.align.roi_file_name];
load2(align_roi_file_name, 'align_roi');

minx = inf;
miny = inf;
maxx = -inf;
maxy = -inf;
tforms = {};
segs_t = {};
xdata = {};
ydata = {};
superpixel_2_seg_map = {};
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('%d ', case_id);
  
  [images, image_prefixes] = get_image_from_stack(config, case_id);

  superpixel_id_offset = 0;
  load2([align_seg_dir, 'sec.', num2str(case_id), '.aligned_seg_mapping', ...
    config.segmentation_choose.choice.seg_suffix, '.mat'], 'label_mappings');
  
  for tile_id = 1:length(image_prefixes)
    if(isempty(image_prefixes{tile_id}))
      continue;
    end
    % load and transform the superpixel segmentation map
    replace_flag = false;
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [superpixel_method, superpixel_suffixes] = ...
      get_superpixel_suffixes(config, image_prefixes{tile_id});
    if(replace_flag)
      config.segmentation_2D = config_segmentation_2D;
    end
   
    superpixel_dir = [get_reconstruction_dir(config), superpixel_config.dir, ...
      superpixel_method, '/'];
    superpixel_seg = load2([superpixel_dir, image_prefixes{tile_id}, ...
      superpixel_suffixes{1}, '.mat']);
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      superpixel_seg.label_map = imresize(superpixel_seg.label_map, ...
        size(images{tile_id}), 'nearest');
    end
    
    load2([sift_dir, image_prefixes{tile_id}, '.sift_transforms.mat']);
    tforms{layer_id, tile_id} = maketform('affine', reshape(transform, [2 3])');
    [segs_t{layer_id, tile_id}, xdata{layer_id, tile_id}, ydata{layer_id, tile_id}] = ...
      imtransform(superpixel_seg.label_map, tforms{layer_id, tile_id}, ...
      'nearest', 'FillValues', -1);
    
    minx = min([minx, xdata{layer_id, tile_id}]);
    miny = min([miny, ydata{layer_id, tile_id}]);
    maxx = max([maxx, xdata{layer_id, tile_id}]);
    maxy = max([maxy, ydata{layer_id, tile_id}]);
    
    % add an offset to the superpixel labels to ensure uniqueness among
    % tiles
    segs_t{layer_id, tile_id}(segs_t{layer_id, tile_id}>0) = ...
      segs_t{layer_id, tile_id}(segs_t{layer_id, tile_id}>0) + superpixel_id_offset;
    
    % for the superpixel_to_seg_label ..
    % load it
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes{tile_id});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    seg_label_map = load2([seg_dir, image_prefixes{tile_id}, ...
      seg_suffix,'.mat']);
    % add superpixel label offset
    seg_label_map.superpixel_to_seg_label(:,1) = ...
      seg_label_map.superpixel_to_seg_label(:,1) + superpixel_id_offset;
    % include mapping of segment ids from the stitching
    seg_label_map.superpixel_to_seg_label(:,2) = ...
      label_mappings{tile_id}(seg_label_map.superpixel_to_seg_label(:,2)+1); %#ok<USENS>
    
    % update the layer's superpixel_2_seg_map
    superpixel_2_seg_map{layer_id}(seg_label_map.superpixel_to_seg_label(:,1)+1) = ...
      seg_label_map.superpixel_to_seg_label(:,2);
    
    superpixel_id_offset = max(segs_t{layer_id, tile_id}(:));
  end
end
fprintf('\n');

for layer_id = 1:size(segs_t,1)
  for tile_id = 1:size(segs_t,2)
    if(isempty(segs_t{layer_id, tile_id}))
      continue;
    end
    xdata{layer_id, tile_id} = round(xdata{layer_id, tile_id}-minx+1);
    ydata{layer_id, tile_id} = round(ydata{layer_id, tile_id}-miny+1);
  end
end
maxx = round(maxx - minx + 1);
maxy = round(maxy - miny + 1);

proofread_config = config.proofreader;
cat = {};
for layer_id = 1:size(segs_t,1)
  case_id = stack_config.case_ids(layer_id);
  fprintf('%d ', case_id);
  
  tile_image = zeros(maxy, maxx)-1;
  for tile_id = 1:size(segs_t,2)
    if(isempty(segs_t{layer_id, tile_id}))
      continue;
    end
    
    temp_image = zeros(maxy, maxx)-1;
    temp_image(ydata{layer_id, tile_id}(1):ydata{layer_id, tile_id}(2), ...
      xdata{layer_id, tile_id}(1):xdata{layer_id, tile_id}(2)) = segs_t{layer_id, tile_id};
    temp_image(tile_image>=0) = tile_image(tile_image>=0);
    tile_image = temp_image;
  end
  tile_image(tile_image<0) = 0;
  tile_image = align_roi.mask .* ...
    tile_image(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  
  if(isfield(proofread_config.roi, 'ymin'))
    tile_image = tile_image(...
      proofread_config.roi.ymin:proofread_config.roi.ymax, ...
      proofread_config.roi.xmin:proofread_config.roi.xmax);
  end
  
  % Remove superpixel ids masked out by roi_mask. Relabel to remove
  % redundant labels.
  superpixel_ids_unique = nonzeros(unique(tile_image(:)));
  superpixel_id_map = zeros(1, max(superpixel_ids_unique)+1);
  superpixel_id_map([1; superpixel_ids_unique+1]) = 0:length(superpixel_ids_unique);
  
  cat{layer_id} = superpixel_id_map(tile_image + 1);
  superpixel_2_seg_map{layer_id} = ...
    [0, superpixel_2_seg_map{layer_id}(superpixel_ids_unique+1)];
  
  % due to upscaling it is possible that some superpixels which were
  % covered during stitching may become exposed. So set any negative
  % superpixel_2_seg_map values to 0
  superpixel_2_seg_map{layer_id}(superpixel_2_seg_map{layer_id}<0) = 0;
end

output_to_proofreader_stitched_superpixel_maps(cat, config);

superpixel_2_seg_map_t = superpixel_2_seg_map; %#ok<NASGU>
save2([get_region_dir(config), 'superpixel_2_seg_map_t', ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], 'superpixel_2_seg_map_t');
fprintf('\ndone.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Merging the linkage graphs from multiple tiles ..\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
links_3D = {};
for layer_id = 1:size(segs_t,1)-1
  case_id = stack_config.case_ids(layer_id);
  case_id_1 = stack_config.case_ids(layer_id+1);
  fprintf('%d ', case_id);
  
  load2([get_reconstruction_dir(config), linkage_config.dir, 'links_3D_p', '.', ...
    num2str(case_id), '.', model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'links_3D_p');
  
  label_mapping_0 = load2([align_seg_dir, 'sec.', num2str(case_id), ...
    '.aligned_seg_mapping', config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'label_mappings');
  label_mapping_1 = load2([align_seg_dir, 'sec.', num2str(case_id_1), ...
    '.aligned_seg_mapping', config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'label_mappings');
  
  links_3D{layer_id} = [];
  for tile_1 = 1:size(links_3D_p, 1)
   for tile_2 = 1:size(links_3D_p,2)
     if(isempty(links_3D_p{tile_1, tile_2}))
       continue;
     end
     links_3D_p{tile_1, tile_2}(:,1) = ...
       label_mapping_0.label_mappings{tile_1}(links_3D_p{tile_1, tile_2}(:,1)+1);
     links_3D_p{tile_1, tile_2}(:,2) = ...
       label_mapping_1.label_mappings{tile_2}(links_3D_p{tile_1, tile_2}(:,2)+1);
     
     links_3D{layer_id} = [links_3D{layer_id}; links_3D_p{tile_1, tile_2}];
   end
  end
  links_3D{layer_id} = unique(sortrows(links_3D{layer_id}), 'rows', 'last');
  
  % Remove segments that have no superpixel members (possibly removed
  % during roi masking)
%   keep_flag = ismember(links_3D_new{layer_id}(:,1), superpixel_2_seg_map{layer_id}(2:end)) ...
%     .* ismember(links_3D_new{layer_id}(:,2), superpixel_2_seg_map{layer_id+1}(2:end));
  keep_flag = links_3D{layer_id}(:,1)>0 & links_3D{layer_id}(:,2)>0;
  links_3D{layer_id} = links_3D{layer_id}(keep_flag,:);
end
save2([get_region_dir(config), 'links_3D', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  config.segmentation_choose.choice.seg_suffix, '.mat'], 'links_3D');

fprintf('\ndone.\n');
return
end
