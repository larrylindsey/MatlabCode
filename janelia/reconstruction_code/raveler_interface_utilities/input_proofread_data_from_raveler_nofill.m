function input_proofread_data_from_raveler_nofill(...
  proofread_data_dir, save_dir, planes, superpixel_map_suffix)
% Get exported proofread data from Raveler

global cat superpixel_2_seg_map proof

superpixel_map_dir = 'sp_maps/';
superpixel_map_prefix = 'sp_map.%d';
if(nargin<=3)
  superpixel_map_suffix = '.tiff';
end
superpixel_2_segment_map_file = 'superpixel-to-segment-map.txt';
segment_2_body_map_file = 'segment-to-body-map.txt';

cat = {};
for p = 1:length(planes)
  plane = planes(p);
  fprintf('plane: %d\n', plane);
  cat{p} = imread(sprintf([proofread_data_dir, superpixel_map_dir, ...
    superpixel_map_prefix, superpixel_map_suffix], plane));
end

superpixel_2_seg_map = {};
fin_sp_2_seg = fopen([proofread_data_dir, superpixel_2_segment_map_file], 'rt');
sp_2_seg = fscanf(fin_sp_2_seg, '%d', [3, inf]);
fclose(fin_sp_2_seg);
for p = 1:length(planes)
  sp_2_seg_p = sp_2_seg(2:3, sp_2_seg(1,:)==p);
  superpixel_2_seg_map{p} = zeros(max(sp_2_seg_p(1,:))+1,1);
  superpixel_2_seg_map{p}(sp_2_seg_p(1,:)+1) = sp_2_seg_p(2,:);
end

proof = [];
fin_seg_2_body = fopen([proofread_data_dir, segment_2_body_map_file], 'rt');
seg_2_body = fscanf(fin_seg_2_body, '%d', [2, inf]);
fclose(fin_seg_2_body);
proof.pmap = zeros(max(seg_2_body(1,:))+1,1);
proof.pmap(seg_2_body(1,:)+1) = seg_2_body(2,:);

%
seg = [];
for i = 1:length(cat)
  seg(:,:,i) = double(proof.pmap( superpixel_2_seg_map{i}(cat{i}+1) + 1)); %#ok<AGROW>
end

%
check_for_dir(save_dir);
save_file_name = [save_dir, 'proofread.seg.mat'];
save(save_file_name, 'seg');
system(['chmod a-w ', save_file_name]);

return;
end

