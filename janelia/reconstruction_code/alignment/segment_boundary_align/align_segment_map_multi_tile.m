function align_segment_map_multi_tile(config)
% align_segment_map_multi_tile(config)
% Align segmentations from multiple tiles within a plane into one
% consistent semgent map.
% The tiles should be aligned before this, for instance, with SIFT-based
% affine transformations.
%
% The segment boundaries are micro-aligned using algorithms such as
% bipartite matching over boundary elements.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

stack_config = config.stack;
align_config = config.align.in_section_align;
align_seg_config = config.align_segmentation;
seg_config = config.segmentation_2D;

if(strcmp(align_config.method, 'SIFT')==1)
  tform_dir = get_sift_dir(config);
  tform_suffix = '.sift_transform.in_plane';
else
  tform_dir = get_deformable_mesh_dir(config);
  tform_suffix = '.dmesh_affine.in_plane';
end
align_seg_dir = get_align_seg_dir(config);

if(~isfield(align_seg_config, 'is_verbose'))
  align_seg_config.is_verbose = true;
end
if(~isfield(align_seg_config, 'is_verbose_figures'))
  align_seg_config.is_verbose_figures = false;
end
if(align_seg_config.is_verbose)
  fprintf('Aligning segment maps of multiple tiles ..\n');
end

for case_id = stack_config.case_ids
  if(align_seg_config.is_verbose)
    fprintf('plane:%d:\n', case_id);
  end
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);

  if(max(is_to_be_processed)==0)
    if(align_seg_config.is_verbose)
      fprintf('All tiles in this section have is_to_be_processed=false, skipping this section\n');
    end
    continue;
  end
  
  % Generate a substring for the file name that accounts for the images
  % involved in the joint alignment
  image_set_string = '.';
  for tile_id = 1:length(image_prefixes)
    image_set_string = [image_set_string, image_prefixes{tile_id}]; %#ok<AGROW>
  end
  image_set_string = strrep(image_set_string, '/', '_');
  image_set_string = strrep(image_set_string, '\', '_');
  if(length(image_prefixes)==1)
    % just one tile in this section so put trivial mapping
    [seg_method, seg_suffix] = ...
      get_segmentation_suffix(config, image_prefixes{1});
    seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
    seg = load2([seg_dir, image_prefixes{1}, seg_suffix, '.mat'], ...
      'label_map');
    label_mappings{1} = 0:max(seg.label_map(:)); %#ok<AGROW>
  else
    % align the segment maps and get segment boundaries
    if(align_seg_config.is_verbose)
      fprintf('Computing transformed segment maps ..\n');
    end
    minx = inf;
    maxx = -inf;
    miny = inf;
    maxy = -inf;
    segs_t = {};
    border_masks = {};
    border_segs = {};
    for tile_id = 1:length(image_prefixes)
      if(align_seg_config.is_verbose)
        fprintf('tile:%d\n', tile_id);
      end
      [seg_method, seg_suffix] = ...
        get_segmentation_suffix(config, image_prefixes{tile_id});
      seg_dir = [get_reconstruction_dir(config), seg_config.dir, seg_method, '/'];
      seg = load2([seg_dir, image_prefixes{tile_id}, seg_suffix, '.mat'], ...
        'label_map');
      fprintf('min(seg.label_map(:)): %d, max(seg.label_map(:)): %d\n', ...
        min(seg.label_map(:)), max(seg.label_map(:)));
      
      load2([tform_dir, image_prefixes{tile_id}, image_set_string, ...
        tform_suffix, '.mat'], 'transform_p');
      if(isfield(config.stack, 'segmentation_scale') && config.stack.segmentation_scale>1)
        transform_p(5:6) = transform_p(5:6)/config.stack.segmentation_scale;
      end
      tform = maketform('affine', reshape(transform_p, [2 3])');

      [segs_t{tile_id}, xdata{tile_id}, ydata{tile_id}] = ...
        imtransform(seg.label_map, tform, 'nearest', 'FillValues', -1); %#ok<AGROW>
      bm = zeros(size(seg.label_map));
      bm(11:end-10,11:end-10) = 1;
      border_masks{tile_id} = ...
        imtransform(bm, tform, 'nearest', 'FillValues', 0); %#ok<AGROW>
      minx = min([minx, xdata{tile_id}]);
      miny = min([miny, ydata{tile_id}]);
      maxx = max([maxx, xdata{tile_id}]);
      maxy = max([maxy, ydata{tile_id}]);

      border_segs{tile_id} = nonzeros(unique([seg.label_map(1,:), seg.label_map(end,:), ...
        seg.label_map(:,1)', seg.label_map(:,end)'])); %#ok<AGROW>
      clear seg
    end
    minx = round(minx);
    miny = round(miny);
    maxx = round(maxx);
    maxy = round(maxy);
    for tile_id = 1:length(image_prefixes)
      if(isempty(image_prefixes{tile_id}))
        continue;
      end
      xdata{tile_id} = round(xdata{tile_id}-minx+1); %#ok<AGROW>
      ydata{tile_id} = round(ydata{tile_id}-miny+1); %#ok<AGROW>
    end

    maxx_0 = maxx - minx + 1;
    maxy_0 = maxy - miny + 1;
    aligned_segs = {};
    aligned_border_masks = {};
    for tile_id = 1:length(image_prefixes)
      if(isempty(image_prefixes{tile_id}))
        continue;
      end
      aligned_segs{tile_id} = zeros(maxy_0, maxx_0)-1; %#ok<AGROW>
      aligned_segs{tile_id}(ydata{tile_id}(1):ydata{tile_id}(2), ...
        xdata{tile_id}(1):xdata{tile_id}(2)) = segs_t{tile_id}; %#ok<AGROW>
      fprintf('min(aligned_segs{tile_id}(:)): %d, max(aligned_segs{tile_id}(:)): %d\n', ...
        min(aligned_segs{tile_id}(:)), max(aligned_segs{tile_id}(:)));
      segs_t{tile_id} = []; %#ok<AGROW>
      if(align_seg_config.is_verbose_figures)
        figure;
        imshow(aligned_segs{tile_id}, []);
        title(['aligned tile segment map: plane ', num2str(case_id), ...
          ' ,tile ', num2str(tile_id)]);
      end
      aligned_border_masks{tile_id} = zeros(maxy_0, maxx_0); %#ok<AGROW>
      aligned_border_masks{tile_id}(ydata{tile_id}(1):ydata{tile_id}(2), ...
        xdata{tile_id}(1):xdata{tile_id}(2)) = border_masks{tile_id}; %#ok<AGROW>
    end
    if(align_seg_config.is_verbose)
      fprintf('done\n');
    end

    % stitch together the segmentations
    if(align_seg_config.is_verbose)
      fprintf('Stitching together segmentations ..\n');
    end
    tile_id = 1;
    if(align_seg_config.is_verbose)
      fprintf('tile:%d\n', tile_id);
    end
    combined_segment_map = aligned_segs{tile_id};
    combined_border_mask = aligned_border_masks{tile_id};
    border_seg = border_segs{tile_id};
    label_mappings = {};
    label_mappings{tile_id} = 0:max(aligned_segs{tile_id}(:)); %#ok<AGROW>
    aligned_segs{tile_id} = []; %#ok<AGROW>
    aligned_border_masks{tile_id} = []; %#ok<AGROW>
    for tile_id = 2:length(image_prefixes)
      if(align_seg_config.is_verbose)
        fprintf('tile:%d\n', tile_id);
      end
      if(isempty(image_prefixes{tile_id}))
        continue;
      end

      switch(align_seg_config.method)
        case 'boundary_bipartite_match'
          segment_correspondence = correspond_segment_boundary_bipartite_match(...
            combined_segment_map, combined_border_mask, border_seg, ...
            aligned_segs{tile_id}, aligned_border_masks{tile_id}, config);
      end

      max_label = 0;
      for t = 1:tile_id-1
        max_label = max(max_label, max(label_mappings{t}));
      end
      fprintf('max_label: %d\n', max_label);
      
      label_offset = max_label + 1;
      aligned_seg_1 = aligned_segs{tile_id};
      aligned_seg_1(aligned_seg_1>0) = aligned_seg_1(aligned_seg_1>0) + label_offset;
      max_label_1 = max(aligned_seg_1(:));
      border_segs{tile_id} = border_segs{tile_id} + label_offset; %#ok<AGROW>
      if(~isempty(segment_correspondence))
        segment_correspondence(2,:) = segment_correspondence(2,:) + label_offset;
      end
      label_mappings{tile_id} = [0, (1:max(aligned_segs{tile_id}(:))) + label_offset]; %#ok<AGROW>
      aligned_segs{tile_id} = []; %#ok<AGROW>
      
      segment_correspondence = [segment_correspondence, ...
        [1:max_label; zeros(1, max_label); zeros(1, max_label)], ...
        [zeros(1, max_label_1); 1:max_label_1; zeros(1, max_label_1)]]; %#ok<AGROW>
      label_mapping = [-1, double(get_pmap_from_linkage_graph({segment_correspondence'}, 1))];

      combined_segment_map_1 = label_mapping(aligned_seg_1+2);
      combined_segment_map_1(combined_segment_map>=0) = ...
        label_mapping(combined_segment_map(combined_segment_map>=0)+2);
      combined_segment_map = combined_segment_map_1;
      clear combined_segment_map_1
      if(align_seg_config.is_verbose_figures)
        random_mapping = randperm(max(combined_segment_map(:)));
        random_mapping = [0, 0, random_mapping]; %#ok<AGROW>
        display_segment_map = random_mapping(combined_segment_map+2);
        figure(3000+tile_id);
        imshow(display_segment_map, []);
      end

%       % find the segments that are completely covered by the previous tiles
%       % and set their mapping to -1
%       max_label = max(label_mapping);
%       label_mapping_2 = [-1 0:max_label];
%       labels_visible = unique(combined_segment_map(:));
%       labels_invisible = ismember(-1:max_label, labels_visible)==0;
%       label_mapping_2(labels_invisible) = -1;
%       label_mapping = label_mapping_2(label_mapping+2);

      % update the segment label mappings
      for t = 1:tile_id
        label_mappings{t} = label_mapping(label_mappings{t}+2); %#ok<AGROW>
      end
      clear label_mapping label_mapping_2 labels_visible labels_invisible
      
      border_seg = [border_seg; border_segs{tile_id}]; %#ok<AGROW>

      combined_border_mask = max(combined_border_mask, aligned_border_masks{tile_id});
      aligned_border_masks{tile_id} = []; %#ok<AGROW>
    end
  end
  clear combined_border_mask combined_segment_map
  
  if(align_seg_config.is_verbose)
    fprintf('done\n');
  end
  
  if(align_seg_config.is_verbose)
    fprintf('Saving ... ');
  end
  save_dir = [align_seg_dir, image_sub_dirs{1}];
  check_for_dir(save_dir);
  save2([save_dir, 'sec.', num2str(case_id), '.aligned_seg_mapping', ...
    image_set_string, config.segmentation_choose.choice.seg_suffix, '.mat'], ...
    'label_mappings');
  if(align_seg_config.is_verbose)
    fprintf('done\n');
  end
end

return
end
