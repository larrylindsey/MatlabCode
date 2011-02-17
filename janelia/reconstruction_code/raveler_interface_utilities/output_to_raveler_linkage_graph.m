function seg_2_body_map = output_to_raveler_linkage_graph(links_3D, superpixel_2_seg_map, ...
  label_pairs_overlap, config)
% output_to_raveler_linkage_graph(links_3D, superpixel_2_seg_map, config)
% Save linkage graph for Raveler.
% If appended_seg_2_body_map is provided, it is appended to the linkage
% graph specified in links_3D.
% 
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%

if(config.is_verbose)
  fprintf('START: output_to_raveler_linkage_graph\n');
end
stack_config = config.stack;
raveler_config = config.proofreader.Raveler;
apply_config = config.linkage.apply;
if(isfield(apply_config, 'linkage_threshold'))
  linkage_threshold = apply_config.linkage_threshold;
else
  linkage_threshold = input('Enter linkage threshold: ');
end

if(isfield(apply_config, 'retain_links_subset'))
  if(config.is_verbose)
    fprintf('Retaining only a subset of links\n');
  end
  switch(apply_config.retain_links_subset)
    case 'max_confidence'
      if(config.is_verbose)
        fprintf('For each segment, link with maximum confidence\n');
      end
      for i = 1:length(links_3D)
        l1 = links_3D{i};
        l1 = sortrows(l1, 3);
        [junk, I] = unique(l1(:,1), 'last'); %#ok<ASGLU>
        l1 = l1(I,:);
        
        l2 = links_3D{i};
        l2 = sortrows(l2, 3);
        [junk, I] = unique(l2(:,2), 'last'); %#ok<ASGLU>
        l2 = l2(I,:);
        
        l3 = links_3D{i}(links_3D{i}(:,3)>2.5, :);
        
        links_3D{i} = [l1; l2; l3];
        links_3D{i} = unique(links_3D{i}, 'rows');
      end
  end
end

% get all segment ids occuring in superpixel_2_seg_map
segment_ids = nonzeros(unique(superpixel_2_seg_map{1}));
links_3D{1} = [links_3D{1}; segment_ids, zeros(size(segment_ids,1),2)];
for i = 2:length(superpixel_2_seg_map)
  seg_ids_l = nonzeros(unique(superpixel_2_seg_map{i}));
  size(seg_ids_l)
  size(links_3D{i-1})
  links_3D{i-1} = [links_3D{i-1}; ...
    zeros(size(seg_ids_l)), seg_ids_l, zeros(size(seg_ids_l))];
  segment_ids = [segment_ids; seg_ids_l]; %#ok<AGROW>
end
segment_ids = unique(segment_ids);


if(isfield(apply_config, 'shrink_wrap') && ...
    apply_config.shrink_wrap.is_enabled ==true)
  fprintf('Shrink-wrap the reconstruction volume\n');
  fprintf('First get the segments touching the perimeter shell\n');
  raveler_config = config.proofreader.Raveler;
  sp_dir = [get_to_be_proofread_dir(config), raveler_config.superpixel_dir];
  case_id = stack_config.case_ids(1);
  sp = imread([sp_dir, 'superpixel_map.', ...
    raveler_config.superpixel_version_name, '.', ...
    num2str(case_id, '%05d'), '.png']);
  perimeter_seg_ids{1} = ...
    uint32(unique(nonzeros(superpixel_2_seg_map{1}(sp+1))));
  for i = 2:length(superpixel_2_seg_map)-1
    case_id = stack_config.case_ids(i);
    fprintf('case_id: %d\n', case_id);
    sp = imread([sp_dir, 'superpixel_map.', ...
      raveler_config.superpixel_version_name, '.', ...
      num2str(case_id, '%05d'), '.png']);
    perimeter_seg_ids{i} = uint32([ ...
      unique(nonzeros(superpixel_2_seg_map{i}(sp(:,1)+1))); ...
      unique(nonzeros(superpixel_2_seg_map{i}(sp(:,end)+1))); ...
      unique(nonzeros(superpixel_2_seg_map{i}(sp(1,:)+1))); ...
      unique(nonzeros(superpixel_2_seg_map{i}(sp(end,:)+1)))]); %#ok<AGROW>
  end
  case_id = stack_config.case_ids(end);
  sp = imread([sp_dir, 'superpixel_map.', ...
    raveler_config.superpixel_version_name, '.', ...
    num2str(case_id, '%05d'), '.png']);
  perimeter_seg_ids{end+1} = ...
    uint32(unique(nonzeros(superpixel_2_seg_map{end}(sp+1))));
  sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy_seeded_seg(...
    links_3D, linkage_threshold, perimeter_seg_ids, ...
    apply_config.shrink_wrap.linkage_threshold);
else
  % Simply use the generated linakge graph as is - no shrink-wrap
  
  % seg_2_body_map = get_pmap_from_linkage_graph_with_dummies(links_3D, ...
  %   linkage_threshold);
%   sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy_sw(...
%     links_3D, label_pairs_overlap, linkage_threshold);
  sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
    links_3D, linkage_threshold);
end

seg_2_body_map = [];
seg_2_body_map(sec_seg_2_body_map(:,2)+1) = sec_seg_2_body_map(:,3);

% make sure all segment ids are being mapped to _some_ body
seg_2_body_map = [seg_2_body_map, length(seg_2_body_map):max(segment_ids)];
  
file_name = [get_to_be_proofread_dir(config), ...
  raveler_config.segment_to_body_map_file_name];
fout_seg_2_body_map = fopen(file_name, 'wt+');
fprintf(fout_seg_2_body_map, '# segment-body map\n# segment	body\n\n');

seg_body_ids = seg_2_body_map(segment_ids+1);
d = int32([0 segment_ids'; 0 seg_body_ids]);
fprintf(fout_seg_2_body_map, '%d\t\t%d\n', d);

fclose(fout_seg_2_body_map);

if(config.is_verbose)
  fprintf('STOP: output_to_raveler_linkage_graph\n');
end
return
end
