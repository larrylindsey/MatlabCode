function body_stack = get_body_label_stack_from_raveler_proofread_data(...
  proofread_data_dir, planes, varargin)
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
superpixel_map_prefix = 'sp_map.%05d';
superpixel_map_suffix = '.png';
superpixel_2_segment_map_file = 'superpixel_to_segment_map.txt';
segment_2_body_map_file = 'segment_to_body_map.txt';
for i = 1:2:length(varargin)
  switch(varargin{i})
    case 'superpixel_map_dir'
      superpixel_map_dir = varargin{i+1};
    case 'superpixel_map_prefix'
      superpixel_map_prefix = varargin{i+1};
    case 'superpixel_map_suffix'
      superpixel_map_suffix = varargin{i+1};
    case 'superpixel_2_segment_map_file'
      superpixel_2_segment_map_file = varargin{i+1};
    case 'segment_2_body_map_file'
      segment_2_body_map_file = varargin{i+1};
  end
end

fin_sp_2_seg = fopen([proofread_data_dir, superpixel_2_segment_map_file], 'rt');
sp_2_seg = fscanf(fin_sp_2_seg, '%d', [3, inf])';
fclose(fin_sp_2_seg);

fin_seg_2_body = fopen([proofread_data_dir, segment_2_body_map_file], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2, inf]);
segment_2_body_map = zeros(max(seg_2_body(1,:))+1,1);
segment_2_body_map(seg_2_body(1,:)+1) = seg_2_body(2,:);
fclose(fin_seg_2_body);

size_sp =   size(imread(sprintf([proofread_data_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], planes(1))));
body_stack = zeros([size_sp, length(planes)]);
for p = 1:length(planes)
  plane = planes(p);
  fprintf('plane: %d\n', plane);

  superpixel_map = imread(sprintf([proofread_data_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], plane));
  
  sp_2_seg_p = sp_2_seg(sp_2_seg(:,1)==plane, 2:3);
  segment_map = apply_mapping(superpixel_map, sp_2_seg_p);
  
  body_map = segment_2_body_map(1+segment_map);
  body_stack(:,:,p) = body_map;
end

body_stack = remove_merged_boundaries_3D(uint32(body_stack));
return;
end
