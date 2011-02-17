
clear all
close all

import_dir = '/groups/chklovskii/chklovskiilab/em_reconstruction/proofread_data/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/region.crop4_global_alignment_0161_1460.unreal.911157326168313_40/ms3_1011.1110_4k.4k_09042009/';
output_tif_stack_file_name = '/groups/visitors/home/takemuras/Desktop/medulla_161-1460/TmX_02946_ds7.tif';
case_ids = 1460:-1:161; %91:-1:1;
xy_scale = 7; %6; %3; %2;
z_scale = 2; %2; %4; % = 6/2

superpixel_map_dir = 'superpixel_maps/';
superpixel_map_prefix = 'sp_map.%05d';
superpixel_map_suffix = '.png';

superpixel_to_segment_file_name = 'superpixel_to_segment_map.txt';
segment_to_body_file_name = 'segment_to_body_map.txt';

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
%  body_id_remap(1+5931) = 255; %L1   4390
%  body_id_remap(1+554) = 255; %L2    1011
%  body_id_remap(1+5559) = 255; %L3   992
%  body_id_remap(1+5954) = 255; %L5   4861
%  body_id_remap(1+5562) = 255; %R7
%  body_id_remap(1+5561) = 255; %R8   1653
%  body_id_remap(1+2992) = 255; %C2
%  body_id_remap(1+2894) = 255; %C3   2218
%  body_id_remap(1+) = 255; %T1
%  body_id_remap(1+2889) = 255; %Tm1   4329
%  body_id_remap(1+2885) = 255; %Tm2   1710 re-examine around #1200-1400
%  body_id_remap(1+2950) = 255; %Tm3   3517
%  body_id_remap(1+2900) = 255; %Mi9   1650
%  body_id_remap(1+2890) = 255; %Tm12? 1804
%  body_id_remap(1+2899) = 255; %Mi4   1647
%  body_id_remap(1+2896) = 255; %Mi1   4056
%  body_id_remap(1+5649) = 255; %Dm8   10284
 body_id_remap(1+2946) = 255; %TmX   5595 neighbor?
%  body_id_remap(1+2895) = 255; %TmX   5489
%  body_id_remap(1+2888) = 255; %TmX   4672
%  body_id_remap(1+440) = 255; %TmX    4672 re-examine
%  body_id_remap(1+2901) = 255; %TmX   4054
%  body_id_remap(1+2928) = 255; %TmX   3965 re-examine
%  body_id_remap(1+2898) = 255; %TmX   3878
%  body_id_remap(1+2892) = 255; %nTm3  5584
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255; %
%  body_id_remap(1+) = 255;
%  body_id_remap(1+) = 255;
%  body_id_remap(1+) = 255;
%  body_id_remap(1+) = 255;
%  body_id_remap(1+) = 255;
%  body_id_remap(1+) = 255;
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
    imwrite(uint8(imresize(uint8(seg_to_body(seg+1)), 1/xy_scale, 'nearest')), ...
      output_tif_stack_file_name, 'WriteMode', 'append', 'Compression', 'none');
  end
end

