
clear all
close all

import_dir = '/groups/visitors/home/takemuras/data/311_610_temp_export_shinya.081709/';
output_tif_stack_file_name = '~/Desktop/medulla_311-610/1582_bodies.tif'; %'~/Desktop/3D viewer/Uo_volume.tif';
case_ids = 610:-1:311; %91:-1:1;
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
%  body_id_remap(1+537) = 1;
%  body_id_remap(1+1516) = 2; %Tm2?
%  body_id_remap(1+1533) = 3;
 body_id_remap(1+1582) = 4;
%  body_id_remap(1+1798) = 5;
%  body_id_remap(1+1893) = 6; %L2 
%  body_id_remap(1+1988) = 7; %Tm
%  body_id_remap(1+3180) = 8; %Tm1
%  body_id_remap(1+3281) = 9;
%  body_id_remap(1+3286) = 10;
%  body_id_remap(1+3290) = 11;
%  body_id_remap(1+3313) = 12;
%  body_id_remap(1+3345) = 13;
%  body_id_remap(1+3631) = 14;
%  body_id_remap(1+3694) = 15;
%  body_id_remap(1+3761) = 16;
%  body_id_remap(1+3835) = 17;
%  body_id_remap(1+3849) = 18;
%  body_id_remap(1+3880) = 19;
%  body_id_remap(1+3882) = 20;
%  body_id_remap(1+3899) = 21;
%  body_id_remap(1+3937) = 22;
%  body_id_remap(1+4033) = 23;
%  body_id_remap(1+4043) = 24;
%  body_id_remap(1+4131) = 25;
%  body_id_remap(1+4189) = 26;
%  body_id_remap(1+4190) = 27;
%  body_id_remap(1+4191) = 28;
%  body_id_remap(1+4195) = 29;
%  body_id_remap(1+4219) = 30;
%  body_id_remap(1+) = 31;
%  body_id_remap(1+) = 32;
%  body_id_remap(1+) = 33;
%  body_id_remap(1+) = 34;
%  body_id_remap(1+) = 35;
%  body_id_remap(1+) = 36;
%  body_id_remap(1+) = 37;
%  body_id_remap(1+) = 38;
%  body_id_remap(1+) = 39;
%  body_id_remap(1+) = 40;
%  body_id_remap(1+) = 41;
%  body_id_remap(1+) = 42; 
%  body_id_remap(1+) = 43;
%  body_id_remap(1+) = 44;
%  body_id_remap(1+) = 45;
%  body_id_remap(1+) = 46;
%  body_id_remap(1+) = 47;
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

