function gen_linkage_gph_intensity_pair_v2_boost_prealign_sift(config)
% gen_linkage_gph_intensity_pair_v2_boost_prealign_sift(config)
% generate linkage graph using boosted classifier trained on
% intensity pair histograms features
% Section alignment :prealigned, single-tile imod, or multi-tile SIFT
% based - see align_stack_SIFT_affine_tile_pair_inter_plane.m.
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
% v2  07152008  appended multiresolution pyramid of the 2D histograms.
%

stack_config = config.stack;
seg_config = config.segmentation_2D(end);
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;

version_names{1} = 'intensity_pair_hist_v2';
version_names{2} = 'intensity_pair_hist_v2b';
version_names{3} = 'intensity_pair_hist_v2c';
version_names{4} = 'intensity_pair_hist_v2d';
version_id = find(strcmp(feature_config.type, version_names));
if(isempty(version_id))
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
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
weak_learner = tree_node_w(linkage_model.tree_depth); %#ok<NASGU>

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n%d ', case_id);

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
label_props_1 = {};
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
  label_props_1{tile} = regionprops(segment_2D_label_map_1{tile}.label_map, 'Area'); %#ok<AGROW>
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  case_id_prev = stack_config.case_ids(i-1);
  fprintf('%d:', case_id);
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
  label_props_2 = {};
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
    label_props_2{tile} = regionprops(segment_2D_label_map_2{tile}.label_map, 'Area'); %#ok<AGROW>
  end

  links_3D_p = {};

  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 1:length(images_1)
    for tile_2 = 1:length(images_2)
      fprintf('%d-%d ', tile_1, tile_2);
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
        catch
          continue;
        end
        if(isempty(transforms_tp))
          continue;
        end
        xdata = {};
        ydata = {};
        tform_1 = maketform('affine', ...
          reshape(transforms_tp(1,:), [2 3])'); %#ok<USENS>
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
          image_1_tt;
        seg_1_tt = imtransform(segment_2D_label_map_1{tile_1}.label_map, tform_1, ...
          'nearest', 'FillValues', -1);
        segment_2D_label_map_1{tile_1}.label_map_t = zeros(maxy, maxx)-1;
        segment_2D_label_map_1{tile_1}.label_map_t...
          (ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2)) = seg_1_tt;

        image_2_t = zeros(maxy, maxx);
        image_2_t(ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2)) = ...
          image_2_tt;
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
            align_tform, margin, 'nearest');

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
          image_1_t = images_1{tile_1};
          image_2_t = images_2{tile_2};
          segment_2D_label_map_1{tile_1}.label_map_t = segment_2D_label_map_1{tile_1}.label_map;
          segment_2D_label_map_2{tile_2}.label_map_t = segment_2D_label_map_2{tile_2}.label_map; %#ok<AGROW>
        end
      end
      
      %%%%%%%%%%%%%%%%%%%%
      % compute intensity pair features
      %%%%%%%%%%%%%%%%%%%%
      % special case for 'intensity_pair_hist_v2c' feature type
      if(version_id==3)
        region_overlap_features = ...
          collect_area_overlap_features_c(segment_2D_label_map_1{tile_1}.label_map_t, ...
          segment_2D_label_map_2{tile_2}.label_map_t);
        
        link_confidence = ...
          Classify(linkage_model.GLearners, linkage_model.GWeights, ...
          region_overlap_features(3:5,:)*stack_config.area_factor);

        links_3D_p{tile_1, tile_2} = [region_overlap_features(1:2,:); link_confidence]'; %#ok<AGROW>
        
        % remove locked segments from the linkage graph
        if(isfield(segment_2D_label_map_1{tile_1}, 'locked_labels'))
          to_retain = ~ismember(links_3D_p{tile_1, tile_2}(:,1), ...
            segment_2D_label_map_1{tile_1}.locked_labels);
          links_3D_p{tile_1, tile_2} = links_3D_p{tile_1, tile_2}(to_retain, :); %#ok<AGROW>
        end
        if(isfield(segment_2D_label_map_2{tile_2}, 'locked_labels'))
          to_retain = ~ismember(links_3D_p{tile_1, tile_2}(:,2), ...
            segment_2D_label_map_2{tile_2}.locked_labels);
          links_3D_p{tile_1, tile_2} = links_3D_p{tile_1, tile_2}(to_retain, :); %#ok<AGROW>
        end
        continue;
      end
      
      % for others
      intensity_pairs = [reshape(image_1_t, [numel(image_1_t), 1]), reshape(image_2_t, ...
        [numel(image_1_t), 1])];
      label_pairs = [reshape(segment_2D_label_map_1{tile_1}.label_map_t, [numel(image_1_t), 1]), ...
        reshape(segment_2D_label_map_2{tile_2}.label_map_t, [numel(image_1_t), 1])];

      intensity_pairs = intensity_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);
      label_pairs = label_pairs(label_pairs(:,1)>0 & label_pairs(:,2)>0, :);

      [unique_label_pairs, I, J] = unique(label_pairs, 'rows');

      %%%
      % build linkage graph for tile pair
      %%%
      links_3D_p{tile_1, tile_2} = []; %#ok<AGROW>
      for p = 1:size(unique_label_pairs, 1)
        intensity_pair_set = intensity_pairs(J==p, :);
        switch(version_id)
          case 1
            h = hist2(intensity_pair_set, linkage_model.intensity_bins, linkage_model.intensity_bins);

            % if the model has been trained of images of resolution different from
            % the current one, then the area and length should be adjusted. This
            % would be especially useful during bootstrapping.
            h = h * stack_config.area_factor;
            h_p = pyramid_hist2(h,2);
            intensity_pair_hist_feature = [reshape(h, [1, numel(h)]), h_p];
          case 2
            h = hist2(intensity_pair_set, linkage_model.intensity_bins, linkage_model.intensity_bins);

            % if the model has been trained of images of resolution different from
            % the current one, then the area and length should be adjusted. This
            % would be especially useful during bootstrapping.
            h = h * stack_config.area_factor;
            h_p = pyramid_hist2(h,2);
            intensity_pair_hist_feature = [label_props_1{tile_1}(unique_label_pairs(p,1)).Area, ...
              label_props_2{tile_2}(unique_label_pairs(p,2)).Area, size(intensity_pair_set, 1), ...
              reshape(h, [1, numel(h)]), h_p];
          case 4
            min_segment_area = min(label_props_1{tile_1}(unique_label_pairs(p,1)).Area, ...
              label_props_2{tile_2}(unique_label_pairs(p,2)).Area);
            intensity_pair_hist_feature = [[label_props_1{tile_1}(unique_label_pairs(p,1)).Area, ...
              label_props_2{tile_2}(unique_label_pairs(p,2)).Area, size(intensity_pair_set, 1), ...
              min_segment_area-size(intensity_pair_set, 1)] * stack_config.area_factor, ...
              size(intensity_pair_set, 1)/min_segment_area];
        end

        features = intensity_pair_hist_feature;
        link_confidence = Classify(linkage_model.GLearners, linkage_model.GWeights, features');

        links_3D_p{tile_1, tile_2} = ...
          [links_3D_p{tile_1, tile_2}; unique_label_pairs(p, :), link_confidence]; %#ok<AGROW>
      end;
    end
  end
  
  save2([get_reconstruction_dir(config), linkage_config.dir, 'links_3D_p', '.', ...
    num2str(case_id_prev), '.', model_config.type, '.', feature_config.type, ...
    apply_config.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
    '.mat'], 'links_3D_p');

  images_1 = images_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  label_props_1 = label_props_2;
  fprintf('\n');
end;

fprintf('done\n');