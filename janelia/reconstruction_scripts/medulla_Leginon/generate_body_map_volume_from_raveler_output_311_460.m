
clear all
close all

import_dir = '/groups/visitors/home/takemurasa/Desktop/Work/311_460_temp_export_satoko.061109/';
output_tif_stack_file_name = '~/Desktop/New_Folder/bundle_12_bodies.tif'; %'~/Desktop/3D viewer/Uo_volume.tif';
case_ids = 460:-1:311; %91:-1:1;
xy_scale = 3; %2;
z_scale = 4; % = 6/2

superpixel_map_dir = 'sp_maps/';
superpixel_map_prefix = 'sp_map.%d';
superpixel_map_suffix = '.png';

superpixel_to_segment_file_name = 'superpixel-to-segment-map.txt';
segment_to_body_file_name = 'segment-to-body-map.txt';

% read in the superpixel-to-segment-map
fin = fopen([import_dir, superpixel_to_segment_file_name], 'rt');
sp_to_seg_map = fscanf(fin, '%d', [3 inf])';
fclose(fin);

% read in the segment-to-body-map
fin = fopen([import_dir, segment_to_body_file_name], 'rt');
seg_to_body_map = fscanf(fin, '%d', [2 inf])';
fclose(fin);
seg_to_body = [];
seg_to_body(1+seg_to_body_map(:,1)) = seg_to_body_map(:,2);

% remap body labels
body_id_remap = zeros(max(seg_to_body)+1, 1);
%  body_id_remap(1+31683) = 1;
%  body_id_remap(1+94458) = 2;
%  body_id_remap(1+3431) = 3;
%  body_id_remap(1+111062) = 4;
%  body_id_remap(1+34992) = 5;
%  body_id_remap(1+3135) = 6;
%  body_id_remap(1+1375) = 7;
%  body_id_remap(1+128110) = 8;
%  body_id_remap(1+98133) = 9;
 body_id_remap(1+20466) = 10;
%  body_id_remap(1+149611) = 11;
%  body_id_remap(1+3215) = 12;
%  body_id_remap(1+112473) = 13;
%  body_id_remap(1+3254) = 14;
%  body_id_remap(1+97977) = 15;
%  body_id_remap(1+101053) = 16;
%  body_id_remap(1+106216) = 17;
%  body_id_remap(1+133566) = 18;
%  body_id_remap(1+6226) = 19;
%  body_id_remap(1+1484) = 20;
%  body_id_remap(1+31986) = 21;
%  body_id_remap(1+104961) = 22;
%  body_id_remap(1+98532) = 23;
%  body_id_remap(1+122690) = 24;
%  body_id_remap(1+95613) = 25;
 body_id_remap(1+5549) = 26;
%  body_id_remap(1+11656) = 27;
%  body_id_remap(1+78260) = 28;
 body_id_remap(1+5360) = 29;
%  body_id_remap(1+17898) = 30;
%  body_id_remap(1+57492) = 31;
 body_id_remap(1+5888) = 32;
%  body_id_remap(1+16826) = 33;
 body_id_remap(1+103809) = 34;
%  body_id_remap(1+2561) = 35;
 body_id_remap(1+98183) = 36;
%  body_id_remap(1+2931) = 37;
%  body_id_remap(1+144360) = 38;
%  body_id_remap(1+114943) = 39;
%  body_id_remap(1+1156774) = 40;
 body_id_remap(1+1128) = 41;
 body_id_remap(1+98145) = 42; 
 body_id_remap(1+1167) = 43;
 body_id_remap(1+115138) = 44;
 body_id_remap(1+1188) = 45;
 body_id_remap(1+1164) = 46;
 body_id_remap(1+6226) = 47;
%  body_id_remap(1+) = 
seg_to_body = body_id_remap(seg_to_body+1);

% read in the superpixel maps
if(exist(output_tif_stack_file_name, 'file')==2)
  delete(output_tif_stack_file_name);
end

for i = 1:length(case_ids)
  case_id = case_ids(i);
  fprintf('case_id: %d\n', case_id);

  sp = imread(sprintf([import_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], case_id));
  sp = imdilate(sp, strel('disk', 2));
  
  sp_to_seg_map_l = sp_to_seg_map(sp_to_seg_map(:,1)==case_id, 2:3);
  sp_to_seg = [];
  sp_to_seg(1+sp_to_seg_map_l(:,1)) = sp_to_seg_map_l(:,2); %#ok<AGROW>
  seg = sp_to_seg(sp+1);

  for zs = 1:z_scale
    imwrite(uint8(imresize(seg_to_body(seg+1), 1/xy_scale, 'nearest')), ...
      output_tif_stack_file_name, 'WriteMode', 'append', 'Compression', 'none');
  end
end

