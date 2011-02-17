function linkage_3D_gen_linkage_gph_intensity_pair_boost(config)
% linkage_3D_gen_linkage_gph_intensity_pair_boost(config)
% generate linkage graph using boosted classifier trained on
% intensity pair histograms features
% 
% - Takes as input a stack of images and 2D superpixel and segment maps.
% - Uses trained classifier from linkage_3D_train_intensity_pair_boost.m to perform 3D linkage.
% - Constructs a linkage graph, links_3D{}, with segments as nodes and 3D link
% confidences as the edge weights. This graph's link-weights are
% thresholded in proofread_linkage_superpixel.m to get the final 3D
% segmentation.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~03202008 init code
% v1  04112008  modified for reconstruction pipeline
%

stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

if(strcmp(feature_config.type, 'intensity_pair_hist')==0)
  error('Feature type does not match with called function. Exiting');
end;
if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

image_dir = get_stack_dir(config);
seg_dir = get_segmentation_dir(config);

junk = tree_node_w(1); %#ok<NASGU>
linkage_model = load([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, '.mat'], ...
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
weak_learner = tree_node_w(linkage_model.tree_depth);

links_3D = [];

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n%d ', case_id);

% load alignment xf data and save for future
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xf_file_name = [image_dir, align_config.xf_file_name];
  xf = read_xf(align_xf_file_name);
  save([image_dir, align_config.xf_save_name], 'xf');
end;

[images_1, image_prefixes_1] = get_image_from_stack(config, case_id);
image_1 = images_1{1};
image_prefix_1 = image_prefixes_1{1};
if(~isempty(stack_config.roi))
  image_1 = image_1(stack_config.roi.ymin:stack_config.roi.ymax, ...
    stack_config.roi.xmin:stack_config.roi.xmax);
end;
image_1 = histeq(im2double(image_1));

segment_2D_label_map_1 = load([seg_dir, image_prefix_1, ...
  apply_config.segmentation_suffix,'.mat']);
if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
  segment_2D_label_map_1.label_map = imresize(segment_2D_label_map_1.label_map, ...
    size(image_1), 'nearest');
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d ', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;
  
  [images_2, image_prefixes_2] = get_image_from_stack(config, case_id);
  image_2 = images_2{1};
  image_prefix_2 = image_prefixes_2{1};
  if(~isempty(stack_config.roi))
    image_2 = image_2(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  image_2 = histeq(im2double(image_2));
  
  segment_2D_label_map_2 = load([seg_dir, image_prefix_2, ...
    apply_config.segmentation_suffix,'.mat']);
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    segment_2D_label_map_2.label_map = imresize(segment_2D_label_map_2.label_map, ...
      size(image_2), 'nearest');
  end
  
  % apply alignment transformation if stack is not pre-aligned
  if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
    align_tform = xf_2_tform(xf(i-1,:), size(image_1, 1), size(image_1, 2));
    margin.left = 1-align_config.margin;
    margin.right = size(image_1,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(image_2,1)+align_config.margin;
    image_1_t = apply_tform(image_1, align_tform, margin, 'bilinear');
    segment_2D_label_map_1.label_map_t = apply_tform(segment_2D_label_map_1.label_map, ...
      align_tform, margin, 'nearest');
    
    align_tform = xf_2_tform(xf(i,:), size(image_2, 1), size(image_2, 2));
    image_2_t = apply_tform(image_2, align_tform, margin, 'bilinear');
    segment_2D_label_map_2.label_map_t = apply_tform(segment_2D_label_map_2.label_map, ...
      align_tform, margin, 'nearest');
    
    if(isfield(apply_config, 'is_verbose') && apply_config.is_verbose)
      figure(898); imshow(image_1_t);
      figure(899); imshow(image_2_t);
      pause;
    end
  else
    image_1_t = image_1;
    image_2_t = image_2;
    segment_2D_label_map_1.label_map_t = segment_2D_label_map_1.label_map;
    segment_2D_label_map_2.label_map_t = segment_2D_label_map_2.label_map;
  end
  
  intensity_pairs = [reshape(image_1_t, [numel(image_1_t), 1]), reshape(image_2_t, ...
    [numel(image_1_t), 1])];
  label_pairs = [reshape(segment_2D_label_map_1.label_map_t, [numel(image_1_t), 1]), ...
    reshape(segment_2D_label_map_2.label_map_t, [numel(image_1_t), 1])];
  
  intensity_pairs = intensity_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
  
  [unique_label_pairs, I, J] = unique(label_pairs, 'rows');
 
  links_3D{i-1} = [];
  for p = 1:size(unique_label_pairs, 1)
    intensity_pair_set = intensity_pairs(J==p, :);
    h = hist2(intensity_pair_set, linkage_model.intensity_bins, linkage_model.intensity_bins);
    
    % if the model has been trained of images of resolution different from
    % the current one, then the area and length should be adjusted. This
    % would be especially useful during bootstrapping.
    h = h * stack_config.area_factor;
    
    intensity_pair_hist_feature = reshape(h, [1, numel(h)]);

    features = [intensity_pair_hist_feature];
    
    link_confidence = Classify(linkage_model.GLearners, linkage_model.GWeights, features');
    
    links_3D{i-1} = [links_3D{i-1}; unique_label_pairs(p, :), link_confidence];
  end;
  
  image_1 = image_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;

%   fprintf('Proceed? (press enter)'); pause;
end;
fprintf('\n');

save([get_reconstruction_dir(config), linkage_config.dir, 'links_3D', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, ...
  apply_config.segmentation_suffix, '.mat'], 'links_3D');

