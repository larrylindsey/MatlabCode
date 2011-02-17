function linkage_evaluate_k_fold_groundtruth_seg_sectionwise(config)
% linkage_evaluate_k_fold_groundtruth_seg_sectionwise(config)
% Train and test on the training data using k-fold protocol - (k-1)/k
% fraction of training data's sections are used for training and rest for
% testing.
%
% The training and testing is performed on the ground-truth segmentation or
% a modified form of it. The labels are assumed to be ground-truth body
% labels (linkage).
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
label_mapping_as_to_gt = {};
for z = 1:n_section
  fprintf('z: %d\n', z);
  imwrite(manual_annotation.al{z}, sprintf([image_dir, stack_config.image_prefix, ...
    stack_config.image_suffix], z));
  
  seg_id = unique(segmentation.wcat{z}(:));
  label_mapping_as_to_gt{z} = []; %#ok<AGROW>
  label_mapping_as_to_gt{z}(1+seg_id) = seg_id; %#ok<AGROW>
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
% total number of segment pairs that were automatically marked to-link
n_to_link = zeros(n_threshold,1);
% total number of segment pairs that were automatically marked to-link that
% are correct
n_to_link_correct = zeros(n_threshold,1);
% total number of segment pairs that were marked to-link in ground-truth
n_to_link_ground_truth = zeros(n_threshold,1);
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
  if(~isfield(evaluate_config, 'is_trained_on_manual_annotation') || ...
      ~evaluate_config.is_trained_on_manual_annotation)
    train_annotation.label_mapping_as_to_gt = {};
  end
  fprintf('Sections included for training:\n');
  for z = [1:(test_bracket_min-1), (test_bracket_max+1):n_section]
    fprintf('%d\n', z);
    train_annotation.al{end+1} = manual_annotation.al{z};
    if(isfield(evaluate_config, 'is_trained_on_manual_annotation') && ...
        evaluate_config.is_trained_on_manual_annotation)
      train_annotation.wcat{end+1} = manual_annotation.wcat{z};
    else
      train_annotation.wcat{end+1} = segmentation.wcat{z};
      train_annotation.label_mapping_as_to_gt{end+1} = ...
        label_mapping_as_to_gt{z};
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
    links_3D_p = linkage_graph.links_3D_p{1,1};
    
    % evaluate linkage confidences
    for t = 1:n_threshold
      linkage_threshold = evaluate_config.link_thresholds(t);

      links_3D_p_t = links_3D_p(links_3D_p(:,3)>=linkage_threshold, :);

      to_link_pairs_al  = unique(links_3D_p_t(:,1:2), 'rows');
      
      n_to_link(t) = n_to_link(t) + size(to_link_pairs_al, 1);
      
      n_to_link_correct(t) = n_to_link_correct(t) + ...
        sum(to_link_pairs_al(:,1)==to_link_pairs_al(:,2));
      
      label_pairs_ground_truth = ...
        unique([manual_annotation.wcat{z}(:), manual_annotation.wcat{z+1}(:)], 'rows');
      label_pairs_ground_truth = label_pairs_ground_truth(...
        min(label_pairs_ground_truth, [], 2)>0, :);
      to_link_pairs_ground_truth = label_pairs_ground_truth(...
        label_pairs_ground_truth(:,1)==label_pairs_ground_truth(:,2), :);
      n_to_link_ground_truth(t) = n_to_link_ground_truth(t) + ...
        size(to_link_pairs_ground_truth, 1);
      
      if(isfield(evaluate_config, 'display_links') && evaluate_config.display_links)
        %%%%%%%%%%%%%%%%
        % Displaying results visually
        %%%%%%%%%%%%%%%%
        
        H = fspecial('gaussian');
        bdry_labels = ones(size(manual_annotation.al{z}));
        bdry_labels(manual_annotation.wcat{z}==0) = 0;
        expanded_bdry_labels = ones(size(bdry_labels));
        expanded_bdry_labels(conv2(bdry_labels,H,'same')<1) = 0;
        centroids1 = regionprops(manual_annotation.wcat{z},'Centroid');
        first_image_wbdry = (255-double(manual_annotation.al{z}))...
          .*expanded_bdry_labels;
        first_image_wbdry(:,:,2) = 255-double(manual_annotation.al{z});
        first_image_wbdry(:,:,3) = (255-double(manual_annotation.al{z}))...
          .*expanded_bdry_labels;
        
        
        
        bdry_labels = ones(size(manual_annotation.al{z}));
        bdry_labels(manual_annotation.wcat{z+1}==0) = 0;
        expanded_bdry_labels = ones(size(bdry_labels));
        expanded_bdry_labels(conv2(bdry_labels,H,'same')<1) = 0;
        
        centroids2 = regionprops(manual_annotation.wcat{z+1},'Centroid');
        
        
        
        second_image_wbdry = (255-double(manual_annotation.al{z+1}))...
          .*expanded_bdry_labels;
        second_image_wbdry(:,:,2) = 255-double(manual_annotation.al{z+1});
        second_image_wbdry(:,:,3) = (255-double(manual_annotation.al{z+1}))...
          .*expanded_bdry_labels;
        
        offset = 255*ones(size(manual_annotation.al{z},1),20);
        
        Big_matrix = [first_image_wbdry(:,:,1),offset,second_image_wbdry(:,:,1)];
        Big_matrix(:,:,2) = [first_image_wbdry(:,:,2),offset,second_image_wbdry(:,:,2)];
        Big_matrix(:,:,3) = [first_image_wbdry(:,:,3),offset,second_image_wbdry(:,:,3)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(1);
        axes('Position', [0 0 1 1]);
        set(gcf, 'Position', [20, 20, size(Big_matrix, 2), size(Big_matrix, 1)]);
        imshow(uint8(Big_matrix));title('Displaying True positive');
        
        true_positive = links_3D_p(links_3D_p(:,3)>= evaluate_config.link_thresholds(t)...
          &  links_3D_p(:,1)==links_3D_p(:,2),1:2);
        
        if size(true_positive,1)>=1
          
          junk = [centroids1(true_positive(1:end,1)).Centroid];
          
          if size(true_positive,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(true_positive,1)]))';
          end
          
          junk = [centroids2(true_positive(1:end,2)).Centroid];
          
          if size(true_positive,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(true_positive,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('True positive Not Found');
        end
        
        saveas(1, [get_reconstruction_dir(config_0), config.linkage.dir, ...
          config.linkage.evaluate.save_dir, 'link_result.tp.', num2str(z), '.', ...
          model_config.type, '.', feature_config.type, ...
          config.linkage.apply.model_suffix, '.tif']);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(2);imshow(uint8(Big_matrix));title('Displaying True Negative');
        
        true_negative = links_3D_p(links_3D_p(:,3)< evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)~=links_3D_p(:,2),1:2);
        
        if size(true_negative,1)>=1
          
          junk = [centroids1(true_negative(1:end,1)).Centroid];
          
          if size(true_negative,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(true_negative,1)]))';
          end
          
          junk = [centroids2(true_negative(1:end,2)).Centroid];
          
          if size(true_negative,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(true_negative,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('True Negative Not Found');
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(3); imshow(uint8(Big_matrix));title('Displaying False positive');
        
        false_positive = links_3D_p(links_3D_p(:,3)>= evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)~=links_3D_p(:,2),1:2);
        
        if size(false_positive,1)>=1
          
          junk = [centroids1(false_positive(1:end,1)).Centroid];
          
          
          if size(false_positive,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(false_positive,1)]))';
          end
          
          junk = [centroids2(false_positive(1:end,2)).Centroid];
          
          if size(false_positive,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(false_positive,1)]))';
          end
          
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('False positive Not Found');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(4); imshow(uint8(Big_matrix));title('Displaying False Negative');
        
        false_negative = links_3D_p(links_3D_p(:,3)< evaluate_config.link_thresholds(t) & ...
          links_3D_p(:,1)==links_3D_p(:,2),1:2);
        
        if size(false_negative,1)>=1
          
          junk = [centroids1(false_negative(1:end,1)).Centroid];
          
          if size(false_negative,1)==1
            st_pt = junk;
          else
            st_pt = (reshape(junk,[2,size(false_negative,1)]))';
          end
          
          junk = [centroids2(false_negative(1:end,2)).Centroid];
          
          if size(false_negative,1)==1
            end_pt = junk;
          else
            end_pt = (reshape(junk,[2,size(false_negative,1)]))';
          end
          end_pt(:,1) = end_pt(:,1)+size(manual_annotation.al{z},2)+20;
          line([(st_pt(:,1))';(end_pt(:,1))'],[(st_pt(:,2))';(end_pt(:,2))'],...
            'linewidth',3);
        else
          title('False negative Not Found');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
    end
  end
  
end

link_eval.n_to_link = n_to_link;
link_eval.n_to_link_correct = n_to_link_correct;
link_eval.n_to_link_ground_truth = n_to_link_ground_truth;

figure;
hPlot = plot(n_to_link_correct./n_to_link_ground_truth, n_to_link_correct./n_to_link, ...
  's-');
hLegend = legend(hPlot, [model_config.type, ' ', feature_config.type, ' ', ...
  config.linkage.apply.model_suffix, ' bp']);
set(hLegend, 'Interpreter', 'none');
title(['Linkage precision-recall ', model_config.type, ' ', ...
  feature_config.type, ' ', config.linkage.apply.model_suffix]);
check_for_dir([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir]);
saveas(gcf, [get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_prec_recall.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.bp.fig']);

save2([get_reconstruction_dir(config_0), config.linkage.dir, ...
  config.linkage.evaluate.save_dir, 'link_prec_recall.', model_config.type, ...
  '.', feature_config.type, config.linkage.apply.model_suffix, '.bp.mat'], 'link_eval');

% save2([get_reconstruction_dir(config), linkage_config.dir, ...
%   'link_dump_k_fold', '.', model_config.type, '.', feature_config.type, ...
%   config.linkage.apply.model_suffix, '.mat'], 'link_dump', '-v7.3');

return
end
