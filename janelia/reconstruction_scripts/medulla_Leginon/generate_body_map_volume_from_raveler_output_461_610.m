
clear all
close all

import_dir = '/groups/visitors/home/takemurasa/Desktop/Work/461_610_temp_export_satoko.062509/';
output_tif_stack_file_name = '~/Desktop/New_Folder/146951_volume.tif'; %'~/Desktop/3D viewer/Uo_volume.tif';
case_ids = 610:-1:461; %91:-1:1;
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
%  body_id_remap(1+6195) = 1;
 body_id_remap(1+146951) = 2;
%  body_id_remap(1+2809) = 3;
%  body_id_remap(1+109944) = 4;
%  body_id_remap(1+71004) = 5;
%  body_id_remap(1+3254) = 6;
%  body_id_remap(1+112514) = 7;
%  body_id_remap(1+122690) = 8;
% body_id_remap(1+112473) = 9;
% body_id_remap(1+98532) = 10;
% body_id_remap(1+2343) = 11;
% body_id_remap(1+12322) = 12;
% body_id_remap(1+18961) = 13;
% body_id_remap(1+5583) = 14;
% body_id_remap(1+12370) = 15;
% body_id_remap(1+5685) = 16;
% body_id_remap(1+6066) = 17;
% body_id_remap(1+6902) = 18;
% body_id_remap(1+5711) = 19;
% body_id_remap(1+29) = 20;
% body_id_remap(1+7112) = 21;
% body_id_remap(1+35279) = 22;
% body_id_remap(1+7002) = 23;
% body_id_remap(1+13344) = 24;
% body_id_remap(1+13428) = 25;
% body_id_remap(1+5549) = 26;
% body_id_remap(1+5600) = 27;
% body_id_remap(1+6309) = 28;
% body_id_remap(1+250) = 29;
% body_id_remap(1+17480) = 30;
% body_id_remap(1+13188) = 31;
% body_id_remap(1+6989) = 32;
% body_id_remap(1+9262) = 33;
% body_id_remap(1+6216) = 34;
% body_id_remap(1+5598) = 35;
% body_id_remap(1+5580) = 37;
% body_id_remap(1+5563) = 38;
% body_id_remap(1+7063) = 39;
% body_id_remap(1+12130) = 36;
%  body_id_remap(1+12647) = 37;
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

