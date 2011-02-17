function linkage_evaluate_k_fold_bodywise(config)
% linkage_evaluate_k_fold(config)
% Train and test on the training data using k-fold protocol - (k-1)/k
% fraction of training data is used for training and rest for testing.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
%

config_0 = config;
stack_config = config.stack;
linkage_config = config.linkage;
feature_config = linkage_config.feature;
model_config = linkage_config.model;
train_config = linkage_config.train;
evaluate_config = linkage_config.evaluate;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(linkage_config.dir, 'dir')~=7)
  mkdir2(linkage_config.dir);
end;
cd(prev_dir);

fprintf('Loading manual annotation .. ');
manual_annotation = load2(train_config.manual_annotation_file);
fprintf('done\n');

n_slice = length(manual_annotation.al);

body_list = [];
for z = 1:n_slice
  body_list = [body_list; nonzeros(unique(manual_annotation.wcat{z}(:)))];
end
body_list = unique(body_list);
n_body = length(body_list);
body_list_orig = body_list;

test_size = round(n_body/evaluate_config.k_fold);

% iteration over the k folds
config.stack.roi = [];
config.stack.align.is_prealigned = true;
config.stack.dir = evaluate_config.temporary_dir;

% save images in a temporary directory
prev_dir = pwd2;
cd(config.stack.dir);
if(exist(stack_config.name, 'dir')~=7)
  mkdir2(stack_config.name);
end;
cd(prev_dir);
image_dir = [config.stack.dir, stack_config.name, '/'];
for z = 1:n_slice
  imwrite(manual_annotation.al{z}, sprintf([image_dir, stack_config.image_prefix, ...
    stack_config.image_suffix], z));
end
config.stack.case_ids = 1:n_slice;
config.superpixel_2_seg.method = 'link_eval_k';
config.segmentation_2D.method = 'link_eval_k';
config.linkage.train.save_suffix = [config_0.linkage.train.save_suffix, 'k'];
config.linkage.apply.model_suffix = config.linkage.train.save_suffix;
config.linkage.apply.segmentation_suffix = evaluate_config.segmentation_suffix;

seg_dir = get_segmentation_dir(config);

% compute linkage confidences, evaluate them and also store in
% link_dump
link_dump = struct('to_be_linked', struct([]), 'not_to_be_linked', struct([]));

n_threshold = length(evaluate_config.link_thresholds);
n_iteration = evaluate_config.n_iteration;
n_total_positive = zeros(n_threshold,n_iteration); % total number of segment pairs that should be linked
n_correct_positive = zeros(n_threshold,n_iteration); % number of segments pairs that were correctly linked
n_total_negative = zeros(n_threshold,n_iteration); % total number of segment pairs that should not be linked
n_correct_negative = zeros(n_threshold,n_iteration); % number of segment pairs that correctly not be linked
for iter = 1:n_iteration
  fprintf('--Iteration %d/%d--\n', iter, n_iteration);
  perm_id = randperm(n_body);
  body_list = body_list_orig(perm_id);
  k = 0;
  for test_bracket_min = 1:test_size:n_body
    k = k+1;
    fprintf('-Fold %d/%d-\n', k, evaluate_config.k_fold);
    test_bracket_max = min(n_body, test_bracket_min+test_size);

    test_body_id = body_list(test_bracket_min:test_bracket_max);
    train_body_id = body_list([1:(test_bracket_min-1), test_bracket_max+1:n_body]);

    %%%%%%%%%%%%%%%
    % Training
    %%%%%%%%%%%%%%%
    % retain only bodies belonging to training set.
    train_annotation = manual_annotation;
    for z = 1:n_slice
      train_annotation.wcat{z} = train_annotation.wcat{z} .* ...
        ismember(train_annotation.wcat{z}, train_body_id);
    end
    save2([image_dir, 'train_annotation_link_eval_k_fold.mat'], '-STRUCT', 'train_annotation');
    pipeline(710, config);

    %%%%%%%%%%%%%%%
    % Testing
    %%%%%%%%%%%%%%%
    % retain only bodies belonging to testing set.
    for z = 1:n_slice
      label_map = manual_annotation.wcat{z} .* ...
        ismember(manual_annotation.wcat{z}, test_body_id); %#ok<NASGU>
      save2(sprintf([seg_dir, stack_config.image_prefix, ...
        evaluate_config.segmentation_suffix,'.mat'] , z), 'label_map');
    end
    pipeline(750, config);

    linkage_graph = load2([get_reconstruction_dir(config), linkage_config.dir, ...
      'links_3D', '.', model_config.type, '.', feature_config.type, ...
      config.linkage.apply.model_suffix, config.linkage.apply.segmentation_suffix, ...
      '.mat'], 'links_3D');
    
    % dump the link pairs for observing the classification
    if(isfield(evaluate_config, 'is_link_dumped') && evaluate_config.is_link_dumped)
      for z = 1:n_slice-1
        to_be_linked_pairs = ...
          find(linkage_graph.links_3D{z}(:,1)==linkage_graph.links_3D{z}(:,2));
        for p = to_be_linked_pairs'
          link_dump.to_be_linked(end+1).z = z;
          link_dump.to_be_linked(end).label_1 = linkage_graph.links_3D{z}(p,1);
          link_dump.to_be_linked(end).label_2 = linkage_graph.links_3D{z}(p,2);
          link_dump.to_be_linked(end).link_conf = linkage_graph.links_3D{z}(p,3);
          [y1, x1] = find(manual_annotation.wcat{z}==linkage_graph.links_3D{z}(p,1));
          [y2, x2] = find(manual_annotation.wcat{z+1}==linkage_graph.links_3D{z}(p,2));
          miny = min([y1; y2]); maxy = max([y1; y2]);
          minx = min([x1; x2]); maxx = max([x1; x2]);
          link_dump.to_be_linked(end).segment_1 = double(manual_annotation.al{z}(miny:maxy, ...
            minx:maxx)) .* double(manual_annotation.wcat{z}(miny:maxy, minx:maxx)== ...
            linkage_graph.links_3D{z}(p,1));
          link_dump.to_be_linked(end).image_1 = double(manual_annotation.al{z}(miny:maxy, ...
            minx:maxx));
          link_dump.to_be_linked(end).segment_2 = double(manual_annotation.al{z+1}(miny:maxy, ...
            minx:maxx)) .* double(manual_annotation.wcat{z+1}(miny:maxy, minx:maxx)== ...
            linkage_graph.links_3D{z}(p,2));
          link_dump.to_be_linked(end).image_2 = double(manual_annotation.al{z+1}(miny:maxy, ...
            minx:maxx));
        end

        not_to_be_linked_pairs = ...
          find(linkage_graph.links_3D{z}(:,1)~=linkage_graph.links_3D{z}(:,2));
        for p = not_to_be_linked_pairs'
          link_dump.not_to_be_linked(end+1).z = z;
          link_dump.not_to_be_linked(end).label_1 = linkage_graph.links_3D{z}(p,1);
          link_dump.not_to_be_linked(end).label_2 = linkage_graph.links_3D{z}(p,2);
          link_dump.not_to_be_linked(end).link_conf = linkage_graph.links_3D{z}(p,3);
          [y1, x1] = find(manual_annotation.wcat{z}==linkage_graph.links_3D{z}(p,1));
          [y2, x2] = find(manual_annotation.wcat{z+1}==linkage_graph.links_3D{z}(p,2));
          miny = min([y1; y2]); maxy = max([y1; y2]);
          minx = min([x1; x2]); maxx = max([x1; x2]);
          link_dump.not_to_be_linked(end).segment_1 = double(manual_annotation.al{z}(miny:maxy, ...
            minx:maxx)) .* double(manual_annotation.wcat{z}(miny:maxy, minx:maxx)== ...
            linkage_graph.links_3D{z}(p,1));
          link_dump.not_to_be_linked(end).image_1 = double(manual_annotation.al{z}(miny:maxy, ...
            minx:maxx));
          link_dump.not_to_be_linked(end).segment_2 = double(manual_annotation.al{z+1}(miny:maxy, ...
            minx:maxx)) .* double(manual_annotation.wcat{z+1}(miny:maxy, minx:maxx)== ...
            linkage_graph.links_3D{z}(p,2));
          link_dump.not_to_be_linked(end).image_2 = double(manual_annotation.al{z+1}(miny:maxy, ...
            minx:maxx));
        end
      end
    end
    
    % evaluate linkage confidences
    for t = 1:n_threshold
      for z = 1:n_slice-1
        n_total_positive(t, iter) = n_total_positive(t, iter) + ...
          sum(linkage_graph.links_3D{z}(:,1)==linkage_graph.links_3D{z}(:,2));
        n_correct_positive(t, iter) = n_correct_positive(t, iter) + ...
          sum(linkage_graph.links_3D{z}(:,3)>= evaluate_config.link_thresholds(t) & ...
          linkage_graph.links_3D{z}(:,1)==linkage_graph.links_3D{z}(:,2));
        n_total_negative(t, iter) = n_total_negative(t, iter) + ...
          sum(linkage_graph.links_3D{z}(:,1)~=linkage_graph.links_3D{z}(:,2));
        n_correct_negative(t, iter) = n_correct_negative(t, iter) + ...
          sum(linkage_graph.links_3D{z}(:,3)< evaluate_config.link_thresholds(t) & ...
          linkage_graph.links_3D{z}(:,1)~=linkage_graph.links_3D{z}(:,2));
      end
    end
  end
end

link_eval.n_total_positive = n_total_positive;
link_eval.n_correct_positive = n_correct_positive;
link_eval.n_total_negative = n_total_negative;
link_eval.n_correct_negative = n_correct_negative;

% figure; plot(n_total_positive-n_correct_positive, n_total_negative-n_correct_negative);

save2([get_reconstruction_dir(config), linkage_config.dir, ...
  'link_eval_k_fold', '.', model_config.type, '.', feature_config.type, ...
  config.linkage.apply.model_suffix, '.mat'], 'link_eval');
save2([get_reconstruction_dir(config), linkage_config.dir, ...
  'link_dump_k_fold', '.', model_config.type, '.', feature_config.type, ...
  config.linkage.apply.model_suffix, '.mat'], 'link_dump', '-v7.3');

return
end
