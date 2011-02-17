function boundary_BEL_dump_testing_data_SAB(config)
% boundary_BEL_dump_testing_data_SAB(config)
% Dump testing data for BEL for running on Windows machines.
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
  fprintf('\nSTART: boundary_BEL_dump_testing_data_SAB\n');
end

bel_dir = [get_reconstruction_dir(config), config.segmentation_2D.dir, ...
  config.BEL.dir];
check_for_dir(bel_dir);

% Initialize job
job_name = ['BEL.test.', num2str(stack_config.case_ids(1)), '.', ...
  num2str(stack_config.case_ids(end))];
job_dir = [get_reconstruction_dir(config), config.job.dir, job_name, '/'];
check_for_dir(job_dir);
check_for_dir([job_dir, 'data/']);
results_dir = ['results_', bel_config.version];
check_for_dir([job_dir, results_dir, '/']);

if(bel_config.dump_test_images)
  fout_backref = fopen([job_dir, job_name, '.image_backref.txt'], 'wt');
  save_dir_name = [job_dir, 'data/'];
  % Dump the images for testing with BEL
  test_image_id = 0;
  for case_id = stack_config.case_ids
    if(bel_config.is_verbose)
      fprintf('plane: %d\n', case_id);
    end
    [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
      get_image_from_stack(config, case_id);

    for tile_id = 1:length(images)
      image = images{tile_id};
      image_prefix = image_prefixes{tile_id};
      if(bel_config.is_verbose)
        fprintf('tile %d: %s\n', tile_id, image_prefix);
      end
      if(is_to_be_processed(tile_id)==0)
        if(bel_config.is_verbose)
          fprintf('is_to_be_processed=false, skipping\n');
        end
        continue;
      end
      
      image = histeq2(image(image>0.1 & image<0.9), image);

      for t = 1:length(bel_config.tiles)
        test_name = [save_dir_name, 'I', num2str(test_image_id, '%05d'), '.tif'];
        fprintf('%s\n', test_name);
        image_crop = image(bel_config.tiles(t).miny:bel_config.tiles(t).maxy, ...
          bel_config.tiles(t).minx:bel_config.tiles(t).maxx);
        imwrite(image_crop, test_name);
        fprintf(fout_backref, '%s %d %d %d %d %s\n', image_prefix, ...
          bel_config.tiles(t).minx, bel_config.tiles(t).maxx, bel_config.tiles(t).miny, ...
          bel_config.tiles(t).maxy, ['data/I', num2str(test_image_id, '%05d'), '.tif']);

        test_image_id = test_image_id + 1;
      end
    end
  end
  fclose(fout_backref);
  n_test_image = test_image_id;
else
  test_image_id = 0;
  for case_id = stack_config.case_ids
    if(bel_config.is_verbose)
      fprintf('plane: %d\n', case_id);
    end
    image_prefixes = get_image_prefixes_subdirs(config, case_id);
    test_image_id = test_image_id + length(image_prefixes)*length(bel_config.tiles);
  end
  n_test_image = test_image_id;
end

% Copy executables for BEL testing
global config_global
copyfile([config_global.bin_dir, 'BEL/*'], [job_dir, '.']);

% Copy BEL training results
copyfile([bel_dir, 'trained_classifiers/', results_dir, '/*'], ...
  [job_dir, results_dir, '/.']);

% Create a BEL config file
BEL_config_file_name = [job_name, '.BEL_config.', bel_config.version, '.txt'];
fout_BEL_config = fopen([job_dir, BEL_config_file_name], 'wt');
fprintf(fout_BEL_config, '%% training or testing\n');
fprintf(fout_BEL_config, '%% 1: training   0: testing\n');
fprintf(fout_BEL_config, 'training=0\n\n');
fprintf(fout_BEL_config, '%% grey sclae or color images\n');
fprintf(fout_BEL_config, '%% 1: color 0: grey scale\n');
fprintf(fout_BEL_config, 'color=%d\n\n', bel_config.color);
fprintf(fout_BEL_config, '%% path to training data\n');
fprintf(fout_BEL_config, 'path source=data\n\n');
fprintf(fout_BEL_config, '%% path to training results\n');
fprintf(fout_BEL_config, 'path results=%s\n\n', results_dir);
fprintf(fout_BEL_config, '%% patch size\n');
fprintf(fout_BEL_config, 'patch size=%d\n\n', bel_config.patch_size);
fprintf(fout_BEL_config, '%% use pre-computed mask, 0: no mask 1: Canny 2: user defined mask\n');
fprintf(fout_BEL_config, 'mask=%d\n\n', bel_config.mask);
fprintf(fout_BEL_config, '%% use Canny detector at different scales as features\n');
fprintf(fout_BEL_config, 'Canny=%d\n\n', bel_config.Canny);
fprintf(fout_BEL_config, 'set2=%d\n\n', bel_config.set2);
fprintf(fout_BEL_config, '%% use position information, e.g. row and col, as candidate features, 0: no, 1: yes\n');
fprintf(fout_BEL_config, 'position=%d\n\n', bel_config.position);
fprintf(fout_BEL_config, '%% to stop training on the node below the portion\n');
fprintf(fout_BEL_config, 'factor to stop=%g\n\n', bel_config.factor_to_stop);
fprintf(fout_BEL_config, '%% number of training/testing images\n');
fprintf(fout_BEL_config, 'number of images=%d\n\n', n_test_image);
fprintf(fout_BEL_config, '%% maximum depth of PBT\n');
fprintf(fout_BEL_config, '%% this largely decides the quality and speed of the final detector\n');
fprintf(fout_BEL_config, 'depth=%d\n\n', bel_config.depth);
fprintf(fout_BEL_config, '%% cascade level\n');
fprintf(fout_BEL_config, 'cascade=%d\n\n', bel_config.cascade);
fprintf(fout_BEL_config, '%% number of weak classifiers in each AdaBoost node\n');
fprintf(fout_BEL_config, 'number of weak classifiers=%d\n\n', ...
  bel_config.number_of_weak_classifiers);
fclose(fout_BEL_config);

% Create a batch script for running BEL testing
fout_job = fopen([job_dir, job_name, '.', bel_config.version,'.bat'], 'wt');
fprintf(fout_job, 'bin\\BEL.exe %s\n', BEL_config_file_name);
fclose(fout_job);

if(bel_config.is_verbose)
  fprintf('\nSTOP: boundary_BEL_dump_testing_data_SAB\n');
end
return
end
