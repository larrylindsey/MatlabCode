
clear all
close all

import_dir = '/groups/chklovskii/chklovskiilab/em_reconstruction/proofread_data/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/region.crop4_global_alignment_0161_0860.161.610.cr79418336731_40/ms3_161.610_3k.3k_combined_stack/';
output_tif_stack_file_name = '/groups/visitors/home/takemuras/Desktop/medulla_161-610/32_bodies.tif';
case_ids = 610:-1:161; %91:-1:1;
xy_scale = 6; %3; %2;
z_scale = 1; %2; %4; % = 6/2

superpixel_map_dir = 'sp_maps/';
superpixel_map_prefix = 'sp_map.%05d';
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
 body_id_remap(1+48) = 128;    %Tm1?
 body_id_remap(1+70) = 129;
 body_id_remap(1+123) = 130;   %postL1
 body_id_remap(1+205) = 131;
 body_id_remap(1+277) = 132;
 body_id_remap(1+324) = 133;
 body_id_remap(1+354) = 134;
 body_id_remap(1+369) = 135;   %R7
 body_id_remap(1+392) = 136;
 body_id_remap(1+483) = 137;
 body_id_remap(1+505) = 138;  %L3
 body_id_remap(1+506) = 139;  %R8
 body_id_remap(1+536) = 140;  %postR8/L3
 body_id_remap(1+730) = 141;  %Tm2
 body_id_remap(1+736) = 142;  %postL3
 body_id_remap(1+738) = 143;
 body_id_remap(1+739) = 144;
 body_id_remap(1+742) = 145;  %postL1
 body_id_remap(1+743) = 146;
 body_id_remap(1+745) = 147;  %postL1/post123
 body_id_remap(1+746) = 148;
 body_id_remap(1+878) = 149;
 body_id_remap(1+920) = 150;
 body_id_remap(1+923) = 151;  %T1?
 body_id_remap(1+925) = 152;  %L1
 body_id_remap(1+930) = 153;  %postL1
 body_id_remap(1+1355) = 154;
 body_id_remap(1+1656) = 155;
 body_id_remap(1+3385) = 156;
 body_id_remap(1+4314) = 157; %L2
 body_id_remap(1+6317) = 158; %postL1
 body_id_remap(1+7808) = 159; %postL1
 body_id_remap(1+29611) = 160;
%  body_id_remap(1+) = 29;
%  body_id_remap(1+) = 30;
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

