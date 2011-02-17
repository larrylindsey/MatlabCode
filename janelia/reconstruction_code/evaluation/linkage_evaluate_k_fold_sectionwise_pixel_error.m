function linkage_evaluate_k_fold_sectionwise_pixel_error(config)
% linkage_evaluate_k_fold_sectionwise_pixel_error(config)
% Train and test on the training data using k-fold protocol - (k-1)/k
% fraction of training data's sections are used for training and rest for
% testing.
%
% Pratim Ghosh,
% Summer intern, Janelia Farm Research Campus, HHMI.
% Univ. of California, Santa Barbara.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus
%
% v1  07092008  init code
% v2  06222009  split data section-wise
% v3  07282009  include modified segmentations.
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

if(isfield(evaluate_config, 'segmentation_file') && ...
    ~isempty(evaluate_config.segmentation_file))
  fprintf('Loading segmentation to be used for training and testing .. ');
  segmentation = load2(evaluate_config.segmentation_file);
  fprintf('done\n');
else
  segmentation = manual_annotation;
end

n_section = length(manual_annotation.al);

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
config.stack.case_ids = 1:n_section;
config.segmentation_2D.method = 'link_eval_k';
config.superpixel_2_seg.method = 'link_eval_k';
config.superpixel_choose.choice.method = 'link_eval_k';
config.superpixel_choose.choice.seg_suffix = '.eval_k';
config.segmentation_choose.choice.method = 'link_eval_k';
config.segmentation_choose.choice.seg_suffix = '.eval_k';
config.linkage.train.save_suffix = [config_0.linkage.train.save_suffix, 'k'];
config.linkage.apply.model_suffix = config.linkage.train.save_suffix;
% config.linkage.apply.segmentation_suffix = evaluate_config.segmentation_suffix;

image_dir = [config.stack.dir, stack_config.name, '/'];
seg_dir = get_segmentation_dir(config);
fprintf('Assigning automatic segmentation to ground-truth segments ...\n');
for z = 1:n_section
  fprintf('z: %d\n', z);
  imwrite(manual_annotation.al{z}, sprintf([image_dir, stack_config.image_prefix, ...
    stack_config.image_suffix], z));
  
  % Assign automatic segmentation map's segments to ground-truth segments.
  % If an A.S. segment overlaps substantially with more than one
  % ground-truth segment then disregard it during training and testing.
  [label_pairs, overlap_area] = count_row_occurence(...
    [segmentation.wcat{z}(:), manual_annotation.wcat{z}(:)]);
  [junk, a_1] = count_row_occurence(segmentation.wcat{z}(:));
  area_1 = [];
  area_1(1+junk) = a_1; %#ok<AGROW>
  [junk, a_2] = count_row_occurence(manual_annotation.wcat{z}(:));
  area_2 = [];
  area_2(1+junk) = a_2; %#ok<AGROW>
  mapping = correspond_labels_many_to_one(label_pairs(:,1), ...
    label_pairs(:,2), overlap_area, area_1', area_2', 200, 0.1);
  label_mapping = zeros(max(segmentation.wcat{z}(:))+1, 1);
  label_mapping(mapping(:,1) + 1) = mapping(:,2);
  
  segmentation.wcat{z} = label_mapping(segmentation.wcat{z}+1);
  
  label_map = segmentation.wcat{z}; %#ok<NASGU>
  save2(sprintf([seg_dir, stack_config.image_prefix, ...
    evaluate_config.segmentation_suffix,'.mat'] , z), 'label_map');
end

config.linkage.train.manual_annotation_file = ...
  [image_dir, 'train_annotation_link_eval_k_fold.mat'];

% compute linkage confidences, evaluate them and also store in
% link_dump
% link_dump = struct('to_be_linked', struct([]), 'not_to_be_linked', struct([]));

n_threshold = length(evaluate_config.link_thresholds);
% total number of pixels linked
n_accept = zeros(n_threshold,1);
% total number of pixels erroneously linked
n_false_accept = zeros(n_threshold,1);
% total number of pixels linked in ground-truth
n_link_gt = zeros(n_threshold,1);
% total number of pixels erroneously left unlinked
n_false_reject = zeros(n_threshold,1);

test_brackets_min = [1, n_section-floor(n_section/2)+1];
test_brackets_max = [floor(n_section/2), n_section];
for k = 1:2
  fprintf('-Fold %d/%d-\n', k, 2);
  test_bracket_min = test_brackets_min(k);
  test_bracket_max = test_brackets_max(k);
  
  %%%%%%%%%%%%%%%
  % Training
  %%%%%%%%%%%%%%%
  fprintf('Training ...\n');
  % retain only bodies belonging to training set.
  train_annotation = [];
  train_annotation.al = {};
  train_annotation.wcat = {};
  fprintf('Sections included for training:\n');
  for z = [1:(test_bracket_min-1), (test_bracket_max+1):n_section]
    fprintf('%d\n', z);
    train_annotation.al{end+1} = manual_annotation.al{z};
    if(isfield(evaluate_config, 'is_trained_on_manual_annotation') && ...
        evaluate_config.is_trained_on_manual_annotation)
      train_annotation.wcat{end+1} = manual_annotation.wcat{z};
    else
      train_annotation.wcat{end+1} = segmentation.wcat{z};
    end
  end
  save2(config.linkage.train.manual_annotation_file, '-STRUCT', 'train_annotation');
  pipeline(710, config);
  fprintf('done.\n');
  
  %%%%%%%%%%%%%%%
  % Generate linkage graph for test sections
  %%%%%%%%%%%%%%%
  fprintf('Generating linkage graphs for test sections ...\n');
  config.stack.case_ids = test_bracket_min:test_bracket_max;
  pipeline(750, config);
  fprintf('done.\n');
  
  %%%%%%%%%%%%%%%
  % Evaluate linkage graph for test sections
  %%%%%%%%%%%%%%%
  fprintf('Evaluating for this fold ...\n');
  for z = test_bracket_min:(test_bracket_max-1)
    fprintf(sprintf('Loading linkage graph for section %d ...\n', z));
    linkage_graph = load2([get_reconstruction_dir(config), linkage_config.dir, ...
      'links_3D_p', '.', num2str(z), '.', model_config.type, '.', feature_config.type, ...
      config.linkage.apply.model_suffix, config.segmentation_choose.choice.seg_suffix, ...
      '.mat'], 'links_3D_p');
    links_3D_p = {linkage_graph.links_3D_p{1,1}};
    
    % make segment ids non-overlaping
    segment_offset = max(links_3D_p{1}(:,1)) + 1;
    links_3D_p{1}(:,2) = links_3D_p{1}(:,2) + segment_offset;
    
    % evaluate linkage confidences
    for t = 1:n_threshold
      linkage_threshold = evaluate_config.link_thresholds(t);
      
      % get the pixelwise body ids
      sec_seg_2_body_map = get_sec_seg_2_body_map_from_links_3D_with_dummy(...
        links_3D_p, linkage_threshold);
      seg_2_body_map = [];
      seg_2_body_map(sec_seg_2_body_map(:,2)+1) = sec_seg_2_body_map(:,3); %#ok<AGROW>
      
      body_map_0 = seg_2_body_map(1 + segmentation.wcat{z});
      wcat_1 = segmentation.wcat{z+1};
      wcat_1(wcat_1>0) = wcat_1(wcat_1>0) + segment_offset;
      body_map_1 = seg_2_body_map(1 + wcat_1);
      
      to_link_map_gt = manual_annotation.wcat{z}==manual_annotation.wcat{z+1} ...
        & manual_annotation.wcat{z}~=0 & manual_annotation.wcat{z+1}~=0;
      to_link_map_as = body_map_0==body_map_1 & body_map_0~=0  & body_map_1~=0;
      
      to_consider_map = manual_annotation.wcat{z}~=0 & ...
        manual_annotation.wcat{z+1}~=0 & body_map_0~=0  & body_map_1~=0;
      false_accept_map = bwareaopen(to_link_map_as & ~to_link_map_gt & ...
        to_consider_map, 100);
      false_reject_map = bwareaopen(~to_link_map_as & to_link_map_gt & ...
        to_consider_map, 100);

      n_accept(t) = n_false_accept(t) + sum(to_link_map_as(:));
      n_false_accept(t) = n_false_accept(t) + sum(false_accept_map(:));
      n_link_gt(t) = n_false_accept(t) + sum(to_link_map_gt(:));
      n_false_reject(t) = n_false_reject(t) + sum(false_reject_map(:));
      
      if(isfield(evaluate_config, 'display_links') && evaluate_config.display_links)
        %%%%%%%%%%%%%%%%
        % Displaying results visually
        %%%%%%%%%%%%%%%%
        display_image_0 = repmat(manual_annotation.al{z}, [1 1 3]);
        display_image_0(:,:,2) = 255*double(to_link_map_gt);
        display_image_1 = repmat(manual_annotation.al{z+1}, [1 1 3]);
        display_image_1(:,:,2) = 255*double(to_link_map_gt);
        display_image = cat(2, display_image_0, ...
          zeros(size(display_image_0,1),10,3), display_image_1);
        figure(1);
        imshow(display_image);
        title('to be linked pixels in ground truth');

        display_image_0 = repmat(manual_annotation.al{z}, [1 1 3]);
        display_image_0(:,:,2) = 255*double(to_link_map_as);
        display_image_1 = repmat(manual_annotation.al{z+1}, [1 1 3]);
        display_image_1(:,:,2) = 255*double(to_link_map_as);
        display_image = cat(2, display_image_0, ...
          zeros(size(display_image_0,1),10,3), display_image_1);
        figure(2);
        imshow(display_image);
        title('to be linked pixels per linkage module');
        
        pause;
      end
    end
  end
  
end

link_eval.n_accept = n_accept;
link_eval.n_link_gt = n_link_gt;
link_eval.n_false_accept = n_false_accept;
link_eval.n_false_reject = n_false_reject;

figure;
hPlot = plot(n_false_reject./n_link_gt, n_false_accept./n_accept, ...
  's-');
axis([0 1 0 1]);
hLegend = legend(hPlot, [model_config.type, ' ', feature_config.type, ' ', ...
  config.linkage.apply.model_suffix]);
set(hLegend, 'Interpreter', 'none');
xlabel('false rejection rate');
ylabel('false acceptance rate');
title(['Linkage fa-fr pixelwise', model_config.type, ' ', ...
  feature_config.type, ' ', config.linkage.apply.model_suffix]);
check_for_dir([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir]);
saveas(gcf, [get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_fa_fr_pixel.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.fig']);

save2([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_eval.fa_fr_pixel.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.mat'], 'link_eval');

return
end
