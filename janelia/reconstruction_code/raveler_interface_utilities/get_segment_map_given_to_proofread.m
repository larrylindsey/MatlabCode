function wcat = get_segment_map_given_to_proofread(...
  data_to_proofread_dir, planes)
% body_stack = get_segment_map_given_to_proofread(...
%   data_to_be_proofread_dir, planes)
% Get the segment maps given for proofreading. This is used in linkage
% evaluation.
% Inputs:
%   data_to_proofread_dir  directory in which data to be proofread is
%                           stored
%   planes                 sequence of planes in the stack
% Output:
%   wcat                  cell array of segmentation maps.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  09082009  init. code
%

superpixel_map_dir = 'superpixel_maps/';
superpixel_map_prefix = 'superpixel_map.v_a.%05d';
superpixel_map_suffix = '.png';
file_name = sprintf([data_to_proofread_dir, superpixel_map_dir, ...
  superpixel_map_prefix, superpixel_map_suffix], planes(1));
if(exist(file_name, 'file')~=2)
  superpixel_map_prefix = 'superpixel_map.%d';
  file_name = sprintf([data_to_proofread_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], planes(1));
  if(exist(file_name, 'file')~=2)
    error('Could not open superpixel map.\n');
  end
end
superpixel_2_segment_map_file = 'superpixel_to_segment_map.txt';

fin_sp_2_seg = fopen([data_to_proofread_dir, superpixel_2_segment_map_file], 'rt');
c = textscan(fin_sp_2_seg, '%d', 'CommentStyle', '#');
fclose(fin_sp_2_seg);
sp_2_seg = reshape(c{1}, [3 length(c{1})/3]);

wcat = cell(1, length(planes));
for p = 1:length(planes)
  plane = planes(p);
  fprintf('plane: %d\n', plane);

  superpixel_map = imread(sprintf([data_to_proofread_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], plane));
  
  sp_2_seg_p = sp_2_seg(2:3, sp_2_seg(1,:)==p);
  superpixel_2_seg_map = zeros(max(sp_2_seg_p(1,:))+1,1);
  superpixel_2_seg_map(sp_2_seg_p(1,:)+1) = sp_2_seg_p(2,:);
  wcat{p} = superpixel_2_seg_map(1+superpixel_map);
end

return;
end
