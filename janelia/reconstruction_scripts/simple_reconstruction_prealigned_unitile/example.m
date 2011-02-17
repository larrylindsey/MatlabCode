function example
%%%
% constants for reading superpixel maps
spmap_dir = '/media/External/Images/Volumejosef-kh_full_alignment/cube_cropped/segment_maps/';
sp_prefix_format = 'segmap_%03d';
sp_suffix = '.png';

boundary_dir = '/media/External/Images/Volumejosef-kh_full_alignment/cube_cropped/prediction/';
boundary_prefix_format = 'prediction_%03d';
boundary_suffix = '.png';
%%%

stack_dir = '/media/External/Images/Volumejosef-kh_full_alignment/cube_cropped/data/';
%stack_dir = '/media/External/Images/Adult R34 Apical Full Volume (final curated 4-1-10)/Volumejosef/Current Images and Traces/';
image_prefix_format = 'image_%03d';
%image_prefix_format = 'R34CA1-B_S12.%d';
%image_suffix = '.jpg';
image_suffix = '.png';
results_dir = '/media/External/Images/Volumejosef-kh_full_alignment/cube_cropped/raveler4/';
%results_dir = '~/temp/3sections_test/';

check_for_dir(results_dir);

section_ids = 45:55;

superpixel_maps = {};
superpixel_to_segment_maps = {};
segment_offset = 0;

minimum_area_threshold = 700;
minimum_boundary_threshold = 0.3;
%minimum_boundary_threshold = 0.1;
%minimum_area_threshold = 2500;


for z = section_ids
  fprintf('z: %d\n', z);
  image_prefix = sprintf(image_prefix_format, z);
  fprintf('image_prefix: %s\n', image_prefix);
  
  image = imread([stack_dir, image_prefix, image_suffix]);
  if size(image, 3) > 1
      image = rgb2gray(image);
  end
  %figure(1);
  %imshow(image);
  %title('image');
  
  fprintf('computing boundary map ...\n');
  
  boundary_image = imread([boundary_dir, ...
      sprintf(boundary_prefix_format, z), boundary_suffix]);
  boundary_map = 1 - im2double(boundary_image);
  boundary_map = medfilt2(boundary_map, [5 5]);
  %figure(2);
  %imshow(boundary_map);
  %title('boundary map');
  
  fprintf('computing superpixel map ...\n');
  seeds = boundary_map<0.3;
  seeds = bwareaopen(seeds, 100, 4);
  boundary_map_s = imimposemin(boundary_map, seeds, 4);
  superpixel_map = watershed(boundary_map_s, 4);
  
  %sp_prefix = sprintf(sp_prefix_format, z);
  %superpixel_map = im2double(imread([spmap_dir, sp_prefix, sp_suffix]));
  
  superpixel_maps{z} = superpixel_map; %#ok<SAGROW>
  %plot_segment_boundaries(image, superpixel_map, 1, 3);
  %title('superpixel map');
  
  fprintf('computing superpixel-to-segment map ...\n');
  [label_map, sp_to_seg] = compute_segmentation_from_superpixels_with_min_area_c(...
    superpixel_map, boundary_map_s, minimum_boundary_threshold, minimum_area_threshold);

  segment_map = double(remove_merged_boundaries_2D(uint32(label_map)));
  
  %sp_to_seg = cat(1, 1:max(superpixel_map), 1:max(superpixel_map))';

  sp_to_seg(:,2) = sp_to_seg(:,2) + segment_offset;
  segment_offset = max(sp_to_seg(:,2));
  
  superpixel_to_segment_maps{z} = sp_to_seg; %#ok<SAGROW>
  %plot_segment_boundaries(image, segment_map, 1, 4);
  %title('segment map');
end

fprintf('performing linkage to RAG ...\n');
links = [];
for z_1 = section_ids(1:end-1)
  z_2 = z_1+1;
  
  image_prefix = sprintf(image_prefix_format, z_1);
  fprintf('image_prefix: %s\n', image_prefix);
  image = imread([stack_dir, image_prefix, image_suffix]);
  if size(image, 3) > 1
      image = rgb2gray(image);
  end
  
  fprintf('computing boundary map 1 ...\n');
  boundary_map_1 = 1 - im2double(image);
  boundary_map_1 = medfilt2(boundary_map_1, [5 5]);

  image_prefix = sprintf(image_prefix_format, z_2);
  fprintf('image_prefix: %s\n', image_prefix);
  image = imread([stack_dir, image_prefix, image_suffix]);
  if size(image, 3) > 1
      image = rgb2gray(image);
  end

  fprintf('computing boundary map 2 ...\n');
  boundary_map_2 = 1 - im2double(image);
  boundary_map_2 = medfilt2(boundary_map_2, [5 5]);
  
  segment_map_1 = apply_mapping(superpixel_maps{z_1}, superpixel_to_segment_maps{z_1});
  segment_map_1 = double(remove_merged_boundaries_2D(uint32(segment_map_1)));
  segment_map_1(boundary_map_1>0.3) = 0;
  segment_map_2 = apply_mapping(superpixel_maps{z_2}, superpixel_to_segment_maps{z_2});
  segment_map_2 = double(remove_merged_boundaries_2D(uint32(segment_map_2)));
  segment_map_2(boundary_map_2>0.3) = 0;
  
  label_pairs_pixels = [segment_map_1(:), segment_map_2(:)];
  [label_pairs, overlap_area] = count_row_occurence(label_pairs_pixels);
  
  % include links
  links = [links; label_pairs, overlap_area]; %#ok<AGROW>
end

links = links(links(:,1)>0 & links(:,2)>0, :);

max_segment_id = segment_offset;

fprintf('building RAG ...\n');
rag = sparse([], [], [], max_segment_id, max_segment_id, 2*size(links, 1));
rag(sub2ind(size(rag), links(:,1), links(:,2))) = links(:,3);
rag(sub2ind(size(rag), links(:,2), links(:,1))) = links(:,3);

area_overlap_threshold = 800;
segment_to_body_map = get_connected_components(rag, area_overlap_threshold);

fprintf('writing images ...\n');
rav_gs_dir = [results_dir, 'grayscale_maps/'];
check_for_dir(rav_gs_dir);
for z = section_ids
  image_prefix = sprintf(image_prefix_format, z);
  fprintf('image_prefix: %s\n', image_prefix);
  image = imread([stack_dir, image_prefix, image_suffix]);
  if size(image, 3) > 1
      image = rgb2gray(image);
  end
  imwrite(image, [rav_gs_dir, 'grayscale.', num2str(z, '%05d'), '.png']);
end

fprintf('writing superpixel map ...\n');
rav_sp_dir = [results_dir, 'superpixel_maps/'];
check_for_dir(rav_sp_dir);
for z = section_ids
  imwrite(uint16(superpixel_maps{z}), [rav_sp_dir, 'sp_map.', num2str(z, '%05d'), '.png'], ...
    'BitDepth', 16);
end

fprintf('writing superpixel to segment map ...\n');
fout = fopen([results_dir, 'superpixel_to_segment_map.txt'], 'wt');
for z = section_ids
  fprintf(fout, '%d\t%d\t%d\n', ...
    [repmat(z, [1 1+ size(superpixel_to_segment_maps{z}, 1)]); ...
    zeros(2, 1), superpixel_to_segment_maps{z}']);
end
fclose(fout);

fprintf('writing segment to body map ...\n');
fout = fopen([results_dir, 'segment_to_body_map.txt'], 'wt');
fprintf(fout, '%d\t%d\n', [zeros(2,1), [1:max_segment_id; segment_to_body_map']]);
fclose(fout);

% your are done
