function body_stack = get_body_label_stack_from_raveler_to_be_proofread_data(...
  proofread_data_dir, planes)
% body_stack = get_body_label_stack_from_raveler_proofread_data(...
%   proofread_data_dir, planes)
% Construct a body label stack from exported proofread data from Raveler.
% Inputs:
%   proofread_data_dir    directory to which data was exported to by
%                           Raveler.
%   planes                sequence of planes in the stack
% Output:
%   body_stack            3D volume of body labels.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  06192009  init. code
%

superpixel_map_dir = 'superpixel_maps/';
superpixel_map_prefix = 'superpixel_map.v1.%05d';
superpixel_map_suffix = '.png'; % '.tiff'
file_name = sprintf([proofread_data_dir, superpixel_map_dir, ...
  superpixel_map_prefix, superpixel_map_suffix], planes(1));
if(exist(file_name, 'file')~=2)
  superpixel_map_prefix = 'superpixel_map.v1.%d';
  file_name = sprintf([proofread_data_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], planes(1));
  if(exist(file_name, 'file')~=2)
    error('Could not open superpixel map.\n');
  end
end
superpixel_2_segment_map_file = 'superpixel_to_segment_map.txt';
segment_2_body_map_file = 'segment_to_body_map.txt';

fin_sp_2_seg = fopen([proofread_data_dir, superpixel_2_segment_map_file], 'rt');
sp_2_seg = fscanf(fin_sp_2_seg, '%d', [3, inf]);
fclose(fin_sp_2_seg);

fin_seg_2_body = fopen([proofread_data_dir, segment_2_body_map_file], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2, inf]);
segment_2_body_map = zeros(max(seg_2_body(1,:))+1,1);
segment_2_body_map(seg_2_body(1,:)+1) = seg_2_body(2,:);
fclose(fin_seg_2_body);

body_stack = [];
for p = 1:length(planes)
  plane = planes(p);
  fprintf('plane: %d\n', plane);

  superpixel_map = imread(sprintf([proofread_data_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], plane));
  
  sp_2_seg_p = sp_2_seg(2:3, sp_2_seg(1,:)==plane);
  superpixel_2_seg_map = zeros(max(sp_2_seg_p(1,:))+1,1);
  superpixel_2_seg_map(sp_2_seg_p(1,:)+1) = sp_2_seg_p(2,:);
  segment_map = superpixel_2_seg_map(1+superpixel_map);
  
  body_map = segment_2_body_map(1+segment_map);
  body_map_d = imdilate(body_map, strel('disk', 3));
  body_map(body_map==0) = body_map_d(body_map==0);
  body_boundary = draw_segment_boundaries_c(double(body_map));
  body_map(body_boundary==1) = 0;
  
  body_stack(:,:,p) = body_map; %#ok<AGROW>
end

return;
end
