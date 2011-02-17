                                        function mitochondria_apply_boosted_detector_intensity_hist(config)
% mitochondria_apply_boosted_detector_intensity_hist(config)
% apply the trained mitochondria detector to confidence maps on image plane
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  ~Feb.2008  init code
% v1  04112008  modified for reconstruction pipeline
%

stack_config = config.stack;
mitochondria_config = config.mitochondria;
train_config = mitochondria_config.train;
feature_config = mitochondria_config.feature;
model_config = mitochondria_config.model;
apply_config = mitochondria_config.apply;

if(strcmp(model_config.type, 'boost')==0)
  error('Classifier type does not match with called function. Exiting');
end;
if(~isfield(model_config, 'version'))
  model_config.version = '';
end

version_names{1} = 'heq_intensity_hist';
version_names{2} = 'intensity_hist';
version_id = find(strcmp(feature_config.type, version_names));
if(isempty(version_id))
  error('Feature type does not match with called function. Exiting');
end;
if(~isfield(feature_config, 'version'))
  feature_config.version = '';
end

load2([get_reconstruction_dir(config), mitochondria_config.dir, train_config.dir, ...
  train_config.model_prefix, '.', model_config.type, model_config.version, '.', ...
  feature_config.type, feature_config.version, '.mat'], 'GLearners', 'GWeights', 'MaxIter', ...
  'window_sizes', 'intensity_bins');

for case_id = stack_config.case_ids
  
  fprintf('case_id: %d\n', case_id);
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);
  
  for tile_id = 1:length(images)
    fprintf('tile_id: %d\n', tile_id);
    img = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
    fprintf('%s\n', image_prefix);
    
    if(is_to_be_processed(tile_id)==0)
      fprintf('is_to_be_processed=false, skipping\n');
      continue;
    end
    
    if(isfield(stack_config, 'roi') && ~isempty(stack_config.roi))
      img = img(stack_config.roi.ymin:stack_config.roi.ymax, ...
        stack_config.roi.xmin:stack_config.roi.xmax);
    end;

    if(version_id==1)
      img = histeq(img);
    end;
    mitochondria_detection_confidence = zeros(size(img))-inf;
    weak_learner = tree_node_w(model_config.tree_depth); %#ok<NASGU>
    for x = 1:500:size(img,2)
      for y = 1:500:size(img,1)
        x1 = min(size(img,2), x+500);
        y1 = min(size(img,1), y+500);
        fprintf('roi(%d:%d,%d:%d)\n', y, y1, x, x1);
        fprintf('collecting samples..\n');
        sample_mask = zeros(size(img));
        sample_mask(y:y1, x:x1) = 1;
        [test_features, sample_loc] = compute_neighbor_intensity_histogram(img, double(sample_mask), ...
          window_sizes, intensity_bins);
        fprintf('done.\n');

        fprintf('classifying samples..\n');
        d = Classify(GLearners, GWeights, test_features');
        fprintf('done.\n');

        sample_loc_id = (sample_loc(:,1)-1)*size(img,1) + sample_loc(:,2);
        mitochondria_detection_confidence(sample_loc_id) = d;
      end
    end
    
    fprintf('saving ..');
    check_for_dir([get_reconstruction_dir(config), mitochondria_config.dir, ...
      image_sub_dir]);
    save_mitochondria_detection(...
      [get_reconstruction_dir(config), mitochondria_config.dir, image_prefix, ...
      apply_config.save_suffix, '.mat'], mitochondria_detection_confidence);
    fprintf('done. ');

    fprintf('\n');
  end
end;

return
end
