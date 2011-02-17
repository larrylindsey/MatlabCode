function gen_linkage_gph_intensity_pair_v2c_boost_deformable_mesh(config)
% gen_linkage_gph_intensity_pair_v2_boost_deformable_mesh(config)
% generate linkage graph using boosted classifier trained on
% intensity pair histograms features
% Section alignment : Deformable mesh - see
% align_stack_deformable_mesh_tile_pair_inter_plane.m
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
% v3  12122008  modified for deformable mesh
%

stack_config = config.stack;
seg_config = config.segmentation_2D(end);
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
apply_config = linkage_config.apply;
dmesh_config = config.align.linkage_align.deformable_mesh;

version_names{1} = 'intensity_pair_hist_v2c';
version_id = find(strcmp(feature_config.type, version_names)); %#ok<EFIND>
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

dmesh_dir = get_deformable_mesh_dir(config);

junk = tree_node_w(1); %#ok<NASGU>
linkage_model = load2([get_reconstruction_dir(config), linkage_config.dir, 'link_model', '.', ...
  model_config.type, '.', feature_config.type, apply_config.model_suffix, '.mat'], ...
  'GLearners', 'GWeights', 'MaxIter', 'tree_depth', 'annotation_file', 'intensity_bins');
weak_learner = tree_node_w(linkage_model.tree_depth); %#ok<NASGU>

case_id = stack_config.case_ids(1);
fprintf('Constructing linkage graph (links_3D) ..\n%d ', case_id);

[images_1, image_prefixes_1, image_sub_dirs_1, is_to_be_processed_1] = ...
  get_image_from_stack(config, case_id);
size_image = size(images_1{1});
segment_2D_label_map_1 = {};
for tile = 1:length(image_prefixes_1)
  [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_1{tile});
  seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
  segment_2D_label_map_1{tile} = load2([seg_dir, image_prefixes_1{tile}, ...
    seg_suffix,'.mat']); %#ok<AGROW>
  if(~isfield(segment_2D_label_map_1{tile}, 'label_map'))
    if(isfield(config, 'segmentation_2D'))
      config_segmentation_2D = config.segmentation_2D;
      replace_flag = true;
    end
    config.segmentation_2D = config.superpixel_2_seg(1);
    [sp_method, sp_suffix] = ...
      get_superpixel_suffixes(config, image_prefixes_1{tile});
    if(replace_flag)
      config.segmentation_2D = config_segmentation_2D;
    end
    sp_suffix = sp_suffix{1};
    sp_dir = [get_reconstruction_dir(config), seg_config.dir, sp_method, '/'];
    superpixel_label_map = load2([sp_dir, image_prefixes_1{tile}, ...
      sp_suffix,'.mat']);
    superpixel_label_map.label_map(superpixel_label_map.label_map<0) = 0;
    sp_to_seg_map = [];
    sp_to_seg_map(segment_2D_label_map_1{tile}.superpixel_to_seg_label(:,1)+1)=...
      segment_2D_label_map_1{tile}.superpixel_to_seg_label(:,2); %#ok<AGROW>
    if(~isempty(sp_to_seg_map))
      segment_2D_label_map_1{tile}.label_map = ...
        sp_to_seg_map(1+superpixel_label_map.label_map); %#ok<AGROW>
    else
      segment_2D_label_map_1{tile}.label_map = []; %#ok<AGROW>
    end
  end
  if(~isempty(segment_2D_label_map_1{tile}.label_map))
    if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
      segment_2D_label_map_1{tile}.label_map = imresize(segment_2D_label_map_1{tile}.label_map, ...
        size_image, 'nearest'); %#ok<AGROW>
    end
    segment_2D_label_map_1{tile}.label_map(segment_2D_label_map_1{tile}.label_map<0) = 0; %#ok<AGROW>
  end
end
% load fold masks
fold_masks_1 = struct('fold_mask', []);
n_patch_1 = zeros(1, length(image_prefixes_1));
fold_dir = get_fold_dir(config);
for tile = 1:length(image_prefixes_1)
  file_name_prefix = [fold_dir, image_prefixes_1{tile}, '.fold_mask'];
  fold_masks_1(tile).fold_mask = load_fold_mask(file_name_prefix, config);
  n_patch_1(tile) = max(fold_masks_1(tile).fold_mask(:));
end

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d:', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;

  [image_prefixes_2, image_sub_dirs_2, is_to_be_processed_2] = ...
    get_image_prefixes_subdirs(config, case_id);
  
  segment_2D_label_map_2 = {};
  for tile = 1:length(image_prefixes_2)
    [seg_method, seg_suffix] = get_segmentation_suffix(config, image_prefixes_2{tile});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    segment_2D_label_map_2{tile} = load2([seg_dir, image_prefixes_2{tile}, ...
      seg_suffix,'.mat']); %#ok<AGROW>
    if(~isfield(segment_2D_label_map_2{tile}, 'label_map'))
      if(isfield(config, 'segmentation_2D'))
        config_segmentation_2D = config.segmentation_2D;
        replace_flag = true;
      end
      config.segmentation_2D = config.superpixel_2_seg(1);
      [sp_method, sp_suffix] = ...
        get_superpixel_suffixes(config, image_prefixes_2{tile});
      if(replace_flag)
        config.segmentation_2D = config_segmentation_2D;
      end
      sp_suffix = sp_suffix{1};
      sp_dir = [get_reconstruction_dir(config), seg_config.dir, sp_method, '/'];
      superpixel_label_map = load2([sp_dir, image_prefixes_2{tile}, ...
        sp_suffix,'.mat']);
      superpixel_label_map.label_map(superpixel_label_map.label_map<0) = 0;
      sp_to_seg_map = [];
      sp_to_seg_map(segment_2D_label_map_2{tile}.superpixel_to_seg_label(:,1)+1)=...
        segment_2D_label_map_2{tile}.superpixel_to_seg_label(:,2); %#ok<AGROW>
      if(~isempty(sp_to_seg_map))
        segment_2D_label_map_2{tile}.label_map = ...
          sp_to_seg_map(1+superpixel_label_map.label_map); %#ok<AGROW>
      else
        segment_2D_label_map_2{tile}.label_map = []; %#ok<AGROW>
      end
    end
    if(~isempty(segment_2D_label_map_2{tile}.label_map))
      if(isfield(stack_config, 'segmentation_scale') && stack_config.segmentation_scale>1)
        segment_2D_label_map_2{tile}.label_map = ...
          imresize(segment_2D_label_map_2{tile}.label_map, ...
          size_image, 'nearest'); %#ok<AGROW>
      end
      segment_2D_label_map_2{tile}.label_map(segment_2D_label_map_2{tile}.label_map<0) = 0; %#ok<AGROW>
    end
  end
  % load fold masks
  fold_masks_2 = struct('fold_mask', []);
  n_patch_2 = zeros(1, length(image_prefixes_2));
  for tile = 1:length(image_prefixes_2)
    file_name_prefix = [fold_dir, image_prefixes_2{tile}, '.fold_mask'];
    fold_masks_2(tile).fold_mask = load_fold_mask(file_name_prefix, config);
    n_patch_2(tile) = max(fold_masks_2(tile).fold_mask(:));
  end

  
  % compute a linkage graph for each pair of overlapping tiles
  for tile_1 = 1:length(image_prefixes_1)
    for tile_2 = 1:length(image_prefixes_2)
      fprintf('--- Tile pair ---\n%d,%d\n', tile_1, tile_2);
      fprintf('%s\n%s\n', image_prefixes_1{tile_1}, image_prefixes_2{tile_2});
      
      if(is_to_be_processed_1(tile_1)==0 && is_to_be_processed_2(tile_2)==0)
        fprintf('Both tiles have is_to_be_processed=false, skipping this pair\n');
        continue;
      end
      
      file_name_prefix = get_file_name_from_tuple(...
        [get_reconstruction_dir(config), linkage_config.dir], ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
      file_name_suffix = ['.', model_config.type, '.', ...
        feature_config.type, apply_config.model_suffix];
      file_name_linkage_graph = [file_name_prefix, file_name_suffix, ...
        config.segmentation_choose.choice.seg_suffix, '.mat'];
      fprintf('Deleting the target linkage file if it already exists\n');
      delete(get_storage_file_name(file_name_linkage_graph));
     
      if(isfield(apply_config, 'delete_linkage_graph_dont_create') && ...
          apply_config.delete_linkage_graph_dont_create)
        continue;
      end

      fprintf('Loading transform from tile 2 to tile 1\n');
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        fprintf('Trying to load transform txt file\n');
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_prefix, '.tforms.txt']);
        transforms_tp = [];
        if(exist(file_name_tforms, 'file')~=2)
          fprintf('Transform txt file does not exist.\n');
        else
          fprintf('Reading transform txt file.\n');
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
            fprintf('Transform txt file was empty.\n');
          else
            fprintf('Successfully read transform txt file\n');
          end
        end
        fprintf('Trying to load transform tif file\n');
        file_name_map = get_storage_file_name([file_name_prefix, '.map.tif']);
        if(exist(file_name_map, 'file')~=2)
          fprintf('Transform tif file does not exist\n');
        else
          transforms_tp.map_mask = uint16(imread(file_name_map));
          fprintf('Successfully read transform tif file\n');
        end
      else
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_2{tile_2}, image_prefixes_1{tile_1}, 'dt.');
        file_name = [file_name_prefix, '.mat'];
        fprintf('Trying to load transform mat file\n');
        transforms_tp = [];
        try
          load2(file_name,'transforms_tp');
          fprintf('Successfully read transform mat file\n');
        catch %#ok<CTCH>
          fprintf('Could not load transform mat file.\n');
        end
      end
      transforms_tp_rev = transforms_tp;

      fprintf('Loading transform from tile 1 to tile 2\n');
      if(isfield(dmesh_config, 'input_tform_format') && ...
          strcmp(dmesh_config.input_tform_format, 'tif_txt'))
        % Read in piecewise transformations as a TIFF file of the map mask
        % and a text file of the affine transformations.
        fprintf('Trying to load transform txt file\n');
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
        file_name_tforms = get_storage_file_name([file_name_prefix, '.tforms.txt']);
        if(exist(file_name_tforms, 'file')~=2)
          fprintf('Transform txt file does not exist.\n');
        else
          fprintf('Reading transform txt file.\n');
          transforms_tp = [];
          fin_tforms = fopen(file_name_tforms, 'rt');
          transforms_tp.transforms = fscanf(fin_tforms, '%g', [6, inf]);
          fclose(fin_tforms);
          if(isempty(transforms_tp.transforms))
            transforms_tp.map_mask = [];
            fprintf('Transform txt file was empty.\n');
          else
            fprintf('Successfully read transform txt file\n');
          end
        end
        fprintf('Trying to load transform tif file\n');
        file_name_map = get_storage_file_name([file_name_prefix, '.map.tif']);
        if(exist(file_name_map, 'file')~=2)
          fprintf('Transform tif file does not exist\n');
        else
          transforms_tp.map_mask = uint16(imread(file_name_map));
          fprintf('Successfully read transform tif file\n');
        end
      else
        % Read in piecewise transformations as a MATLAB .mat file
        file_name_prefix = get_file_name_from_tuple(dmesh_dir, ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'dt.');
        file_name = [file_name_prefix, '.mat'];
        fprintf('Trying to load transform mat file\n');
        try
          load2(file_name,'transforms_tp');
          fprintf('Successfully read transform mat file\n');
        catch %#ok<CTCH>
          fprintf('Could not load transform mat file. Going to next pair\n');
        end
      end
      
      links_3D_p = [];
      
      if(~isfield(apply_config, 'compute_new_linkage_graph') || ...
        apply_config.compute_new_linkage_graph)
        region_overlap_features = [];
        region_overlap_features_rev = [];
        if(~( isempty(segment_2D_label_map_1{tile_1}.label_map) || ...
            isempty(segment_2D_label_map_2{tile_2}.label_map) || ...
            isempty(transforms_tp) || ...
            isempty(transforms_tp.transforms) || ...
            isempty(transforms_tp.map_mask) || ...
            size(transforms_tp.transforms, 2)<=4))
          
          fprintf('Computing linkage features between segment maps ...\n');
          region_overlap_features = ...
            collect_area_overlap_features_piecewise_affine_map_c(...
            segment_2D_label_map_1{tile_1}.label_map, ...
            segment_2D_label_map_2{tile_2}.label_map, ...
            transforms_tp.map_mask, transforms_tp.transforms);
        end
        if(~( isempty(segment_2D_label_map_1{tile_1}.label_map) || ...
            isempty(segment_2D_label_map_2{tile_2}.label_map) || ...
            isempty(transforms_tp_rev) || ...
            isempty(transforms_tp_rev.transforms) || ...
            isempty(transforms_tp_rev.map_mask) || ...
            size(transforms_tp_rev.transforms, 2)<=4))
          
          fprintf('Computing linkage features between segment maps in reverse...\n');
          region_overlap_features_rev = ...
            collect_area_overlap_features_piecewise_affine_map_c(...
            segment_2D_label_map_2{tile_2}.label_map, ...
            segment_2D_label_map_1{tile_1}.label_map, ...
            transforms_tp_rev.map_mask, transforms_tp_rev.transforms);
        end
        
        if(~isempty(region_overlap_features_rev))
          region_overlap_features = [region_overlap_features, ...
            region_overlap_features_rev([2 1 4 3 5], :)]; %#ok<AGROW>
        end
        
        if(~isempty(region_overlap_features) && ...
            ~(size(region_overlap_features,2)==1 && ...
            max(region_overlap_features(:))==0))
          fprintf('Computing linkage confidences\n');
          link_confidence = ...
            Classify(linkage_model.GLearners, linkage_model.GWeights, ...
            region_overlap_features(3:5,:)*stack_config.area_factor);
          
          links_temp = [region_overlap_features(1:2,:); link_confidence]';
          links_temp = sortrows(links_temp);
          links_temp = unique(links_temp, 'rows');
          
          % remove locked segments from the linkage graph
          if(isfield(segment_2D_label_map_1{tile_1}, 'locked_labels'))
            to_retain = ~ismember(links_temp(:,1), ...
              segment_2D_label_map_1{tile_1}.locked_labels);
            links_temp = links_temp(to_retain, :);
          end
          if(isfield(segment_2D_label_map_2{tile_2}, 'locked_labels'))
            to_retain = ~ismember(links_temp(:,2), ...
              segment_2D_label_map_2{tile_2}.locked_labels);
            links_temp = links_temp(to_retain, :);
          end
          
          links_3D_p = links_temp;
          fprintf('number of overlap pairs: %d\n', size(links_3D_p,1));
        end
      end
      
      % Include proofread linkage graph if asked to do so.
      if(isfield(linkage_config, 'include_proofread') && ...
          linkage_config.include_proofread.is_enabled)
        fprintf('Including proofread linkage graph, proofread_name: %s\n', ...
          linkage_config.include_proofread.proofread_name);
        file_name_prefix = get_file_name_from_tuple(...
          [get_reconstruction_dir(config), linkage_config.dir], ...
          image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'lkp.');
        file_name = [file_name_prefix, '.proofread', '.', ...
          linkage_config.include_proofread.proofread_name, ...
          linkage_config.include_proofread.model_suffix, ...
          config.segmentation_choose.choice.seg_suffix, '.mat'];
        linkage_graph_proofread = [];
        try
          fprintf('Attempting to load %s ...\n', file_name); 
          linkage_graph_proofread = load2(file_name, 'links_3D_p');
        catch %#ok<CTCH>
          fprintf('failed.\nProceeding nevertheless.\n');
        end
        if(~isempty(linkage_graph_proofread))
          fprintf('successfully loaded.\n');
          % Remove all references to proofread segments in links_3D_p and
          % add linkage_graph_proofread to it.
          if(~isempty(links_3D_p))
            to_retain = ~ismember(links_3D_p(:,1), ...
              linkage_graph_proofread.links_3D_p(:,1));
            links_3D_p = links_3D_p(to_retain, :);
          end
          if(~isempty(links_3D_p))
            to_retain = ~ismember(links_3D_p(:,2), ...
              linkage_graph_proofread.links_3D_p(:,2));
            links_3D_p = links_3D_p(to_retain, :);
          end
          fprintf('Added proofread linkage graph to links_3D_p.\n');
          links_3D_p = [links_3D_p; linkage_graph_proofread.links_3D_p]; %#ok<AGROW,NASGU>
        end
      else
        fprintf('However, linkage graph was empty - skipping.\n');
      end
      
      fprintf('saving linkage graph:\n%s\n', ...
        get_storage_file_name(file_name_linkage_graph));  
      save2(file_name_linkage_graph, 'links_3D_p');
    end
  end
  
  image_prefixes_1 = image_prefixes_2;
  segment_2D_label_map_1 = segment_2D_label_map_2;
  is_to_be_processed_1 = is_to_be_processed_2;
  fprintf('\n');
end;

fprintf('done\n');
