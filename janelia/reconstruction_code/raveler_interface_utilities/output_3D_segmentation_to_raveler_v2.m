function output_3D_segmentation_to_raveler_v2(image_tif, image_filter_param, ws_raw, sp_raw, seg_raw, output_dir)
% output_3D_segmentation_to_raveler_v2(image_tif, image_filter, ws_raw, ...
%   sp_raw, seg_raw, output_dir)
% Output the result of 3D segmentation to Raveler such that the superpixels
% and segments are in 2D.
%
% Inputs:
%   image_tif     3D tif stack of the images. Give empty ('') if grayscale
%                   images need not be exported.
%   image_filter  String specifying any filtering to be done on the image.
%                   If none required then give empty ('').
%   ws_raw        Raw file of the 3D watershed
%   sp_raw        Raw file of the mapping from watershed to superpixels
%   seg_raw       Raw file of the mapping from watershed to segments. It is
%                   assumed that the segments are strictly groupings of the
%                   superpixels. If not then an error is returned.
%   output_dir    Directory to which the reconstruction is to be exported
%                   for proofreading.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%


fprintf('START: output_3D_segmentation_to_raveler\n');

check_for_dir(output_dir);

[ws, err] = fread_raw_array_mex(ws_raw);
if(err~=0)
  error('Could not read watershed file [%d]', err);
end

if(~isempty(image_tif))
  fprintf('Dumping grayscale maps ...\n');
  grayscale_dir = [output_dir, 'grayscale_maps/'];
  system(['rm -rf ', grayscale_dir, '*.png']);
  check_for_dir(grayscale_dir);
  for i = 1:size(ws,3)
    fprintf('i: %d\n', i);
    image = im2double(imread(image_tif, i));
    if(size(image,1)~=size(ws,2) || size(image,2)~=size(ws,1))
      error('Dimensions of image and watershed don''t match');
    end
    if(~isempty(image_filter_param))
      image = filter_image2(image, image_filter_param);
    end
    imwrite(image, sprintf([grayscale_dir, 'image.v1.%05d.png'], i));
  end
end

fprintf('Reading in 3D superpixel map\n');
[a, err] = fread_raw_array_mex(sp_raw);
if(err~=0)
  error('Error reading superpixel mapping file [%d]', err);
end
a = a';
m = [];
m(1+a(:,1)) = a(:,2);
superpixel_map = uint32(m(1+ws));
superpixel_map = remove_merged_boundaries_3D(superpixel_map);

fprintf('Reading in 3D segment map\n');
[a, err] = fread_raw_array_mex(seg_raw);
if(err~=0)
  error('Error reading segmentation mapping file [%d]', err);
end
a = a';
m = [];
m(1+a(:,1)) = a(:,2);
segment_map = uint32(m(1+ws));
segment_map = remove_merged_boundaries_3D(segment_map);

fprintf('Dumping superpixel maps ...\n');
superpixel_2_segment_map = [];
segment_2_body_map = [0 0];
superpixel_dir = [output_dir, 'superpixel_maps/'];
system(['rm -rf ', superpixel_dir, '*.png']);
check_for_dir(superpixel_dir);
segment_offset = 0;
for i = 1:size(ws,3)
  fprintf('i: %d\n', i);
  sp_z_3D = superpixel_map(:,:, i)';
  seg_z_3D = segment_map(:,:, i)';
  
  sp_z_2D = relabel_connected_components_2D(sp_z_3D);

  seg_z_2D = relabel_connected_components_2D(seg_z_3D);
  
  body_z_3D = seg_z_3D;
  
  seg_z_2D(seg_z_2D>0) = seg_z_2D(seg_z_2D>0) + segment_offset;
  segment_offset = max(segment_offset, max(seg_z_2D(:)));
  
  sp_2_seg_z_2D = unique([sp_z_2D(:), seg_z_2D(:)], 'rows');
  sp_2_seg_z_2D = sp_2_seg_z_2D(min(sp_2_seg_z_2D, [], 2)>0, :);
  sp_2_seg_z_2D = [0 0; sp_2_seg_z_2D]; %#ok<AGROW>
  
  seg_2_body_z_2D = unique([seg_z_2D(:), body_z_3D(:)], 'rows');
  seg_2_body_z_2D = seg_2_body_z_2D(min(seg_2_body_z_2D, [], 2)>0, :);
  seg_2_body_z_2D = [0 0; seg_2_body_z_2D]; %#ok<AGROW>
  
  [~,f] = unique(sp_2_seg_z_2D(:,1), 'first');
  [~,l] = unique(sp_2_seg_z_2D(:,1), 'last');
  if(nnz(f~=l)>0)
    error('Segments are not strict agglomerations of superpixels');
  end
  [~,f] = unique(seg_2_body_z_2D(:,1), 'first');
  [~,l] = unique(seg_2_body_z_2D(:,1), 'last');
  if(nnz(f~=l)>0)
    error('Bodies are not strict agglomerations of segments');
  end
  
  imwrite(uint16(sp_z_2D), ...
    sprintf([superpixel_dir, 'superpixel_map.v1.%05d.png'], i), ...
    'BitDepth', 16);
  
  superpixel_2_segment_map = [superpixel_2_segment_map;
    repmat(i, [size(sp_2_seg_z_2D,1), 1]), sp_2_seg_z_2D]; %#ok<AGROW>
  
  segment_2_body_map = [segment_2_body_map; seg_2_body_z_2D]; %#ok<AGROW>
end

segment_2_body_map = unique(segment_2_body_map, 'rows');

fprintf('Dumping superpixel to segment map\n');
fout = fopen([output_dir, 'superpixel_to_segment_map.txt'], 'wt');
fprintf(fout, '%d %d %d\n', superpixel_2_segment_map');
fclose(fout);

fprintf('Dumping segment to body map\n');
fout = fopen([output_dir, 'segment_to_body_map.txt'], 'wt');
fprintf(fout, '%d %d\n', segment_2_body_map');
fclose(fout);

fprintf('STOP: output_3D_segmentation_to_raveler\n');
return
end
