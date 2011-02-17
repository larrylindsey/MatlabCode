function gen_linkage_gph_overlap_area_boost_prealign_sift(config)
% gen_linkage_gph_overlap_area_boost_prealign_sift(config)
% generate linkage graph using boosted classifier trained on
% area of overlap features
% Section alignment :prealigned, single-tile imod, or multi-tile SIFT
% based - see align_stack_SIFT_affine_tile_pair_inter_plane.m.
%
% - Takes as input a stack of images and 2D superpixel and segment maps.
% - Uses trained classifier from linkage_3D_train_overlap_area_boost.m to
% perform 3D linkage. 
% - Constructs a linkage graph, links_3D{}, with segments as nodes and 3D link
% confidences as the edge weights. This graph's link-weights are
% thresholded in proofread_linkage_superpixel.m to get the final 3D
% segmentation.
%
% Pratim Ghosh,
% Summer intern, Janelia Farm Research Campus, HHMI.
% Univ. of California, Santa Barbara.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  06202009 init code
%

fprintf('START: gen_linkage_gph_overlap_area_boost_prealign_sift\n');
stack_config = config.stack;
seg_config = config.segmentation_2D(end);
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

if(strcmp(feature_config.type, 'overlap_area')==0)
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
sift_dir = get_sift_dir(config);

junk = tree_node_w(1); %#ok<NASGU>
linkage_model = load2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, '.mat'], ...
  'classifier', 'n_dimension', 'model_config', 'annotation_file');
weak_learner = tree_node_w(linkage_model.model_config.tree_depth); %#ok<NASGU>

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n');

% load alignment xf data and save for future
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xf_file_name = [image_dir, align_config.xf_file_name];
  xf = read_xf(align_xf_file_name);
  save2([image_dir, align_config.xf_save_name], 'xf');
end;

[images_1, image_prefixes_1] = get_image_from_stack(config, case_id);
if(~isempty(stack_config.roi))
  for tile = 1:length(images_1)
    images_1{tile} = images_1{tile}(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end
end;
segment_2D_label_map_1 = {};
for tile = 1:length(images_1)
  images_1{tile} = histeq(im2double(images_1{tile}));
  [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_1{tile});
  seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
  segment_2D_label_map_1{tile} = load2([seg_dir, image_prefixes_1{tile}, ...
    seg_suffix,'.mat']); %#ok<AGROW>
  if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
    segment_2D_label_map_1{tile}.label_map = imresize(segment_2D_label_map_1{tile}.label_map, ...
      size(images_1{tile}), 'nearest'); %#ok<AGROW>
  end
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  case_id_prev = stack_config.case_ids(i-1);
  fprintf('plane_1: %d, plane_2: %d\n', case_id_prev, case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;

  [images_2, image_prefixes_2] = get_image_from_stack(config, case_id);
  if(~isempty(stack_config.roi))
    for tile = 1:length(images_2)
      images_2{tile} = images_2{tile}(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end
  end;
  segment_2D_label_map_2 = {};
  for tile = 1:length(images_2)
    images_2{tile} = histeq(im2double(images_2{tile}));
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_2{tile});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    segment_2D_label_map_2{tile} = load2([seg_dir, image_prefixes_2{tile}, ...
      seg_suffix,'.mat']); %#ok<AGROW>
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      segment_2D_label_map_2{tile}.label_map = ...
        imresize(segment_2D_label_map_2{tile}.label_map, ...
        size(images_2{tile}), 'nearest'); %#ok<AGROW>
    end
  end

  links_3D_p = {};

  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 1:length(images_1)
    for tile_2 = 1:length(images_2)
      fprintf('tile_1: %d, tile_2: %d\n', tile_1, tile_2);
      %%%%%%%%%%%%%%%%%%%%
      % For aligning the two images and segmentation maps
      %%%%%%%%%%%%%%%%%%%%
      if(isfield(stack_config, 'image_structure') && ...
          ~isempty(stack_config.image_structure) ) % multi tile system
        file_name_suffix = get_file_name_from_tuple(sift_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'st.');
        file_name = [file_name_suffix, '.mat'];
        transforms_tp = [];
        try
          load2(file_name,'transforms_tp');
        catch %#ok<CTCH>
          continue;
        end
        if(isempty(transforms_tp))
          continue;
        end
        xdata = {};
        ydata = {};
        tform_1 = maketform('affine', ...
          reshape(transforms_tp(1,:), [2 3])');
        [image_1_tt, xdata{1}, ydata{1}] = ...
          imtransform(images_1{tile_1}, tform_1); %#ok<AGROW>
        minx = min(xdata{1}); miny = min(ydata{1});
        maxx = max(xdata{1}); maxy = max(ydata{1});

        tform_2 = maketform('affine', ...
          reshape(transforms_tp(2,:), [2 3])');
        [image_2_tt, xdata{2}, ydata{2}] = ...
          imtransform(images_2{tile_2}, tform_2); %#ok<AGROW>
        minx = min([minx, xdata{2}]);
        miny = min([miny, ydata{2}]);
        maxx = max([maxx, xdata{2}]);
        maxy = max([maxy, ydata{2}]);
        for t = 1:2
          xdata{t} = round(xdata{t}-minx+1); %#ok<AGROW>
          ydata{t} = round(ydata{t}-miny+1); %#ok<AGROW>
        end
        maxx = round(maxx - minx + 1);
        maxy = round(maxy - miny + 1);

        image_1_t = zeros(maxy, maxx);
        image_1_t(ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2)) = ...
          image_1_tt; %#ok<NASGU>
        seg_1_tt = imtransform(segment_2D_label_map_1{tile_1}.label_map, tform_1, ...
          'nearest', 'FillValues', -1);
        segment_2D_label_map_1{tile_1}.label_map_t = zeros(maxy, maxx)-1; %#ok<AGROW>
        segment_2D_label_map_1{tile_1}.label_map_t...
          (ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2)) = seg_1_tt; %#ok<AGROW>

        image_2_t = zeros(maxy, maxx);
        image_2_t(ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2)) = ...
          image_2_tt; %#ok<NASGU>
        seg_2_tt = imtransform(segment_2D_label_map_2{tile_2}.label_map, tform_2, ...
          'nearest', 'FillValues', -1);
        segment_2D_label_map_2{tile_2}.label_map_t = zeros(maxy, maxx)-1; %#ok<AGROW>
        segment_2D_label_map_2{tile_2}.label_map_t...
          (ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2)) = seg_2_tt; %#ok<AGROW>
      else
        % apply alignment transformation if stack is not pre-aligned
        if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
          align_tform = xf_2_tform(xf(i-1,:), size(images_1{tile_1}, 1), ...
            size(images_1{tile_1}, 2));
          margin.left = 1-align_config.margin;
          margin.right = size(images_1{tile_1},2)+align_config.margin;
          margin.top = 1-align_config.margin;
          margin.bottom = size(images_1{tile_1},1)+align_config.margin;

          image_1_t = apply_tform(images_1{tile_1}, align_tform, margin, 'bilinear');
          segment_2D_label_map_1{tile_1}.label_map_t = ...
            apply_tform(segment_2D_label_map_1{tile_1}.label_map, ...
            align_tform, margin, 'nearest'); %#ok<AGROW>

          align_tform = xf_2_tform(xf(i,:), size(images_2{tile_2}, 1), ...
            size(images_2{tile_2}, 2));
          image_2_t = apply_tform(images_2{tile_2}, align_tform, margin, 'bilinear');
          segment_2D_label_map_2{tile_2}.label_map_t = ...
            apply_tform(segment_2D_label_map_2{tile_2}.label_map, ...
            align_tform, margin, 'nearest'); %#ok<AGROW>

          if(isfield(apply_config, 'is_verbose') && apply_config.is_verbose)
            figure(898); imshow(image_1_t);
            figure(899); imshow(image_2_t);
            pause;
          end
        else
          image_1_t = images_1{tile_1}; %#ok<NASGU>
          image_2_t = images_2{tile_2}; %#ok<NASGU>
          segment_2D_label_map_1{tile_1}.label_map_t = segment_2D_label_map_1{tile_1}.label_map; %#ok<AGROW>
          segment_2D_label_map_2{tile_2}.label_map_t = segment_2D_label_map_2{tile_2}.label_map; %#ok<AGROW>
        end
      end
      
      seg_0 = segment_2D_label_map_1{tile_1}.label_map_t;
      seg_1 = segment_2D_label_map_2{tile_2}.label_map_t;
      label_pairs = [...
        reshape(seg_0, [numel(seg_0), 1]), reshape(seg_1, [numel(seg_1), 1])];
      label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
      [linkage_label_pairs_unique, linkage_label_pairs_count] = ...
        count_row_occurence(label_pairs);
      [s_1, a_1] = count_row_occurence(seg_0(:));
      areas_1 = [];
      areas_1(s_1+1) = a_1; %#ok<AGROW>
      [s_2, a_2] = count_row_occurence(seg_1(:));
      areas_2 = [];
      areas_2(s_2+1) = a_2; %#ok<AGROW>
      features = [areas_1(1+linkage_label_pairs_unique(:,1))', ...
        areas_2(1+linkage_label_pairs_unique(:,2))', linkage_label_pairs_count];
      
      if(size(features,2)~=linkage_model.n_dimension)
        error('Feature dimensions during training and testing don''t match.');
      end
      link_confidence = Classify(linkage_model.classifier.GLearners, ...
        linkage_model.classifier.GWeights, features'*stack_config.area_factor);
      
      links_3D_p{tile_1, tile_2} = ...
        [linkage_label_pairs_unique, link_confidence']; %#ok<AGROW>
    end
  end
  
  save2([get_reconstruction_dir(config), linkage_config.dir, 'links_3D_p', '.', ...
    num2str(case_id_prev), '.', model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'links_3D_p');

  images_1 = images_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  fprintf('\n');
end;

fprintf('done\n');

fprintf('STOP: gen_linkage_gph_overlap_area_boost_prealign_sift\n');
return
end
