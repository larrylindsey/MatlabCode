function boundary_BEL_copy_test_results_to_boundary_maps(config)
% boundary_BEL_copy_test_results_to_boundary_maps(config)
% Copy testing results from BEL to boundary maps
%
% Input:
%   config    config datastructure of the reconstruction
%
% Boosted Edge Learner binary (BEL.exe) was provided by Piotr Dollar et al.
% Please see code/bin/BEL/ReadMe.txt for lisencing information. The program
% is based on Piotr Dollar, Zhuowen Tu, and Serge Belongie, "Supervised
% Learning of Edges and Object Boundaries", CVPR, June, 2006.
% This program is purely for non-profit research purposes. The files within
% code/bin/BEL should not be distributed without the consent of Piotr
% Dollar, Zhuowen Tu, and Serge Belongie.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI.
%
% v0  03132009  init. code
%

stack_config = config.stack;
bel_config = config.precompute.boundary.BEL;

if(~isfield(bel_config, 'is_verbose'))
  dmesh_config.is_verbose = true;
end
if(bel_config.is_verbose)
  fprintf('\nSTART: boundary_BEL_copy_test_results_to_boundary_maps\n');
end

bel_dir = [get_reconstruction_dir(config), config.segmentation_2D.dir, ...
  config.BEL.dir];

job_name = ['BEL.test.', num2str(stack_config.case_ids(1)), '.', ...
  num2str(stack_config.case_ids(end))];
job_dir = [get_reconstruction_dir(config), config.job.dir, job_name, '/'];
results_dir = ['results_', bel_config.version];

% read in the back reference file
back_ref.image_prefixes = {};
back_ref.tiles = [];
back_ref.bel_test_prefixes = {};
fin_backref = fopen([job_dir, job_name, '.image_backref.txt'], 'rt');
s = fgets(fin_backref);
i = 0;
while(s~=-1)
  i = i + 1;
  b = strsplit(s, ' ');
  back_ref.image_prefixes{i} = b{1};
  back_ref.tiles(i,:) = [str2double(b{2}), str2double(b{3}), str2double(b{4}), ...
    str2double(b{5})];
  back_ref.bel_test_prefixes{i} = get_file_name_from_full_path(b{6});
  s = fgets(fin_backref);
end
fclose(fin_backref);

% Get the results from BEL, combine into boundary maps, and save them.
images = get_image_from_stack(config, stack_config.case_ids(1));
[height, width] = size(images{1});
for case_id = stack_config.case_ids
  if(bel_config.is_verbose)
    fprintf('plane: %d\n', case_id);
  end
  [image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_prefixes_subdirs(config, case_id);

  for tile_id = 1:length(image_prefixes)
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
    if(bel_config.is_verbose)
      fprintf('tile %d: %s\n', tile_id, image_prefix);
    end
    if(is_to_be_processed(tile_id)==0)
      if(bel_config.is_verbose)
        fprintf('is_to_be_processed=false, skipping\n');
      end
      continue;
    end
    
    test_id = find(strcmp(image_prefix, back_ref.image_prefixes));
    if(isempty(test_id))
      error('Can not find BEL results for image!\n');
    end
    boundary = ones(height, width);
    for t = test_id
      bel_result_name = [job_dir, results_dir, '/', ...
        back_ref.bel_test_prefixes{t}(1:end-4), '_prob.tif'];
      fprintf('%s\n', bel_result_name);
      bel_result = im2double(imread(bel_result_name));
      bel_result = filter_image(bel_result, bel_config.post_process_filter);
      temp = ones(height, width);
      temp(back_ref.tiles(t,3):back_ref.tiles(t,4), ...
        back_ref.tiles(t,1):back_ref.tiles(t,2)) = bel_result;
      boundary = min(boundary, temp);
    end
    check_for_dir([bel_dir, image_sub_dir]);
    save2([bel_dir, image_prefix, '.BEL.', bel_config.version, ...
      '_', bel_config.post_process_filter, '.mat'], 'boundary');
  end
end

if(bel_config.is_verbose)
  fprintf('\nSTOP: boundary_BEL_copy_test_results_to_boundary_maps\n');
end
return
end
