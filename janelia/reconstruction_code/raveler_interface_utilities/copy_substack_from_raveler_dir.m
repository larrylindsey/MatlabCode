function copy_substack_from_raveler_dir(source, target, case_ids)

check_for_dir(target);

fprintf('copying grayscale maps\n');
o_dir = [target, 'grayscale_maps/'];
check_for_dir(o_dir);
system(['rm ', o_dir, '*']);
for i = 1:length(case_ids)
  case_id = case_ids(i);
  fprintf('case_id: %d\n', case_id);
  comm = ['cp ', source, 'grayscale_maps/', ...
    sprintf('image.v_a.%05d.png', case_id), ' ', o_dir, '.'];
  system(comm);
end

fprintf('copying superpixel maps\n');
o_dir = [target, 'superpixel_maps/'];
check_for_dir(o_dir);
system(['rm ', o_dir, '*']);
for i = 1:length(case_ids)
  case_id = case_ids(i);
  fprintf('case_id: %d\n', case_id);
  comm = ['cp ', source, 'superpixel_maps/', ...
    sprintf('superpixel_map.v_a.%05d.png', case_id), ' ', o_dir, '.'];
  system(comm);
end

fprintf('creating superpixel to segment map\n');
comm = ['rm ', target, 'superpixel_to_segment_map.txt'];
system(comm);
for i = 1:length(case_ids)
  case_id = case_ids(i);
  fprintf('case_id: %d\n', case_id);
  comm = ['grep "^', num2str(case_id), '[^0-9]" ', ...
    source, 'superpixel_to_segment_map.txt | uniq >> ', ...
    target, 'superpixel_to_segment_map.txt'];
  system(comm);
end

fprintf('creating segment to body map\n');
if(true)
  comm = ['rm ', target, 'segment_to_body_map.txt'];
  system(comm);
  fprintf('Sorting superpixel to segment map\n');
  comm = ['more +4 ', source, 'superpixel_to_segment_map.txt  | ', ...
    'sort -n +2 -3  > ', target, 'temp1.txt'];
  system(comm);
  fprintf('Sorting segment to body map\n');
  comm = ['more +4 ', source, 'segment_to_body_map.txt | ', ...
    'sort -n +0 -1  > ' target, 'temp2.txt'];
  system(comm);
  fprintf('Joining superpixel-to-segment and segment-to-body maps\n');
  comm = ['join -1 3 -2 1 ', target, 'temp1.txt ', target, 'temp2.txt', ...
    ' 2> /dev/null | awk ''{print $2 "\t" $3 "\t" $1 "\t" $4;}'' > ', target, 'temp3.txt'];
  system(comm);
  fprintf('Retaining required segment mappings\n');
  for i = 1:length(case_ids)
    case_id = case_ids(i);
    fprintf('case_id: %d\n', case_id);
    comm = ['grep "^', num2str(case_id), '[^0-9]" ', ...
      target, 'temp3.txt | awk ''{print $3 "\t" $4;}'' | uniq >> ', ...
      target, 'segment_to_body_map.txt'];
    system(comm);
  end
else
  comm = ['grep -E ''"''`cat ', target, 'superpixel_to_segment_map.txt | awk ''{print $3;}'' | uniq | awk ''{print "^"$1"[^0-9]";}'' | tr ''\n'' "|"`foogoo''"'' ', ...
    source, 'segment_to_body_map.txt > ', target, 'segment_to_body_map.txt'];
  system(comm);
end
return
end
