function mitochondria_3D_linkage_confidence(config)
% mitochondria_3D_linkage_confidence(config)
% Compute linkage wieghts between segments in adjacent layers by computing
% their percent area overlap after erosion by 3.
%
% Nicholas Sofroniew,
% Univ. of Cambridge, UK.
% Visitor March 2009, JFRC, HHMI.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  03122009  init code
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
segment_config = mitochondria_config.segment;
classify_config = segment_config.classify;

if(~isfield(mitochondria_config, 'is_verbose'))
  mitochondria_config.is_verbose = true;
end

mito_dir = [get_reconstruction_dir(config), mitochondria_config.dir];

if(mitochondria_config.is_verbose)
  fprintf('\nSTART: mitochondria_3D_linkage_confidence\n');
end

for i = 1:length(stack_config.case_ids)
  case_id_0 = stack_config.case_ids(i);
  if(mitochondria_config.is_verbose)
    fprintf('plane:%d\n', case_id_0);
  end
  if(mitochondria_config.is_verbose)
    fprintf('Loading segments ... ');
  end
  image_prefixes = get_image_prefixes_subdirs(config, case_id_0);
  image_prefix_0 = image_prefixes{1};

  load_file_name = [mito_dir, image_prefix_0, ...
    '.mito.seg', segment_config.mito_seg_suffix,...
    '.mat'];
  segments_0 = load2(load_file_name, 'label_map');

  if(i<length(stack_config.case_ids))
    case_id_u = stack_config.case_ids(i+1);
    image_prefixes = get_image_prefixes_subdirs(config, case_id_u);
    image_prefix_u = image_prefixes{1};
    load_file_name = [mito_dir, image_prefix_u, ...
      '.mito.seg', segment_config.mito_seg_suffix,...
      '.mat'];
    segments_u = load2(load_file_name, 'label_map');
  end
  if(i>1)
    case_id_d = stack_config.case_ids(i-1);
    image_prefixes = get_image_prefixes_subdirs(config, case_id_d);
    image_prefix_d = image_prefixes{1};
    load_file_name = [mito_dir, image_prefix_d, ...
      '.mito.seg', segment_config.mito_seg_suffix,...
      '.mat'];
    segments_d = load2(load_file_name, 'label_map');
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  StEl=strel('disk',3); %% Size of erosion applied before computing
  % percent area overlap
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%% TURN THIS INTO A FILTER THAT CAN BE APPLIED %%%%%%%%%%%

  if(mitochondria_config.is_verbose)
    fprintf('Performing erosions ... ');
  end
  eroded_segs_0=imerode(segments_0.label_map,StEl);
  if(i<length(stack_config.case_ids))
    eroded_segs_u=imerode(segments_u.label_map,StEl);
  end
  if(i>1)
    eroded_segs_d=imerode(segments_d.label_map,StEl);
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  if(mitochondria_config.is_verbose)
    fprintf('Computing linkage features ... ');
  end
  num_segs=max(segments_0.label_map(:));
  linkage_features=[];
  %%% For each segment compute its linkage features
  for seg_index=1:num_segs
    region=eroded_segs_0==seg_index;
    area=sum(sum(region));
    if (area~=0)
      if(i<length(stack_config.case_ids))
        intersect=unique(eroded_segs_u(region));
        for intesect_index=2:length(intersect)
          area_overlap=length(find(eroded_segs_u(region)==intersect(intesect_index)));
          area_link_seg=length(find(eroded_segs_u==intersect(intesect_index)));
          linkage_features=[linkage_features; seg_index,intersect(intesect_index),1, ...
            area,area_link_seg,area_overlap/area];
        end
      end
      if(i>1)
        intersect=unique(eroded_segs_d(region));
        for intesect_index=2:length(intersect)
          area_overlap=length(find(eroded_segs_d(region)==intersect(intesect_index)));
          area_link_seg=length(find(eroded_segs_d==intersect(intesect_index)));
          linkage_features=[linkage_features; seg_index,intersect(intesect_index),0,...
            area,area_link_seg,area_overlap/area];
        end
      end
    end
  end
  if(mitochondria_config.is_verbose)
    fprintf('done.\n');
  end

  %% Write matrix with segment features for each image
  save_file_name = [mito_dir, image_prefix_0, ...
    '.mito.seg', segment_config.mito_seg_suffix, ...
    '.linkage_conf', ...
    '.mat'];
  save2(save_file_name, 'linkage_features');

end;
if(mitochondria_config.is_verbose)
  fprintf('\nSTOP: mitochondria_3D_linkage_confidence\n');
end
return;
end
