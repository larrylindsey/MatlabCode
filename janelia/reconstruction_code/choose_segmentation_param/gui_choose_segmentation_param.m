function superpixel_param_choices = gui_choose_segmentation_param(config)
% superpixel_param_choices = gui_choose_segmentation_param(config)
% GUI for choosing segmentation parameters from a user-defined set
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  111908  init code: user specified binary decision tree
% v1  112008  show two segmentations on same image with different colors
% v2  112008  show two segmentations on an image and flip between the two;
%               no original image.
%

stack_config = config.stack;
seg_config = config.segmentation_2D;
choose_config = config.segmentation_choose;
tree_config = choose_config.decision_tree;

if(tree_config.type~=2)
  error('Error: Tree type not binary');
end

superpixel_param_choice_dir = get_superpixel_param_choice_dir(config);

seg_dir = [get_reconstruction_dir(config), seg_config.dir];

superpixel_param_choices = cell(1,length(stack_config.case_ids));
index = 1;
tile_id = 1;
while(index<=length(stack_config.case_ids))
  case_id = stack_config.case_ids(index);
  fprintf('%d:\n', case_id);
  
  [images, image_prefixes, image_sub_dirs, is_to_be_processed] = ...
    get_image_from_stack(config, case_id);
  
  if(tile_id==-1)
    tile_id = length(images);
  end
  while(tile_id<=length(images))
    fprintf('tile %d ', tile_id);
    if(is_to_be_processed(tile_id)==0)
      fprintf('is_to_be_processed=false, skipping\n');
      tile_id = tile_id + 1;
      continue;
    end
    
    image = images{tile_id};
    if(isfield(choose_config, 'is_image_histeq') && choose_config.is_image_histeq)
      image = histeq(image);
    end;
    image_prefix = image_prefixes{tile_id};
    image_sub_dir = image_sub_dirs{tile_id};
    
    if(isfield(stack_config, 'segmentation_scale') && ...
        ~isempty(stack_config.segmentation_scale))
      image = imresize(image, 1/stack_config.segmentation_scale, ...
        stack_config.segmentation_scaling_method);
    end
    if(choose_config.display_top_half)
      image = image(1:round(size(image,1)/2),:);
    end
    
    current_node_id = 1;
    while(current_node_id>0)
      seg = struct('label_map', []);
      param_id = zeros(1,2);
      for i = 1:2
        param_id(i) = tree_config.node(current_node_id).param_id(i);
        seg(i) = load2([seg_dir, choose_config.param(param_id(i)).method, '/', ...
          image_prefix, choose_config.param(param_id(i)).seg_suffix, '.mat'], 'label_map');
        if(choose_config.display_top_half)
          seg(i).label_map = seg(i).label_map(1:round(size(seg(i).label_map,1)/2),:);
        end
      end
      h_fig = figure(1);
      clf;
      set(h_fig, 'Name', ...
        ['Image, plane: ', num2str(case_id), ' tile ', num2str(tile_id), ...
        '.   Choice: ', num2str(param_id(1)), ' and ', num2str(param_id(2))]);
      axes('Position', [0 0 1 1]);
      set(gcf, 'Position', [20, 20, 2600, 1500]);
      user_choice = -1;
      while(user_choice<=0)
        imshow(image);
        colormap('gray');
        hold on;
        h_text = text(round(size(image,2)/2)-250, -20, ...
          ['Choice: ', num2str(param_id(1)), ' and ', num2str(param_id(2)), '. ', ...
          'Press "c" to flip views, "', ...
          choose_config.input_keys{4}, '" to choose ', num2str(param_id(1)), ...
          ', "', choose_config.input_keys{5}, '" to choose ', num2str(param_id(2)), ...
          ', "', choose_config.input_keys{1}, '" to abort', ...
          ', "', choose_config.input_keys{2}, '" redo previous image' ...
          ', "', choose_config.input_keys{3}, '" to redo current image']);
        set(h_text, 'FontSize', 20);
        if(user_choice==-1)
          for i = 1:2
            [py, px] = find(seg(i).label_map==0);
            plot(px, py, ['.', choose_config.plot_colors{i}], 'MarkerSize', 5);
          end
        elseif(user_choice==-2)
          [py, px] = find(seg(2).label_map==0);
          plot(px, py, ['.', choose_config.plot_colors{2}], 'MarkerSize', 5);
        end
        hold off;
        press_type = waitforbuttonpress;
        while(press_type==0)
          press_type = waitforbuttonpress;
        end
        key_press = get(h_fig, 'CurrentCharacter');

        if(strcmp(key_press, 'c')==1)
          % for flipping between views
          if(user_choice==-1)
            user_choice = -2;
          else
            user_choice = -1;
          end
        else
          % choose a parameter
          for k = 1:2+3
            if(strcmp(key_press, choose_config.input_keys{k})==1)
              user_choice = k;
              break;
            end
          end
        end
      end
      switch(user_choice)
        case 1
          error('Program aborted');
        case 2
          % redo previous image
          if(tile_id==1)
            % goto last tile of previous section
            if(index>1)
              index = index - 1;
              tile_id = -1;
              current_node_id = 1;
              break;
            end
          else
            % goto previous tile in current section
            tile_id = tile_id - 1;
            current_node_id = 2;
            break;
          end
        case 3
          % redo current image
          current_node_id = 1;
        otherwise
          % progress down decision tree
          current_node_id = tree_config.node(current_node_id).node_id(user_choice-3);
      end
    end

    % For current_node_id, leaves of decision tree have negative node ids,
    % whereas positive values indicate exceptions.
    choice_param_id = -current_node_id;
    
    % For choice_param_id, positive value indicates a valid segmentation
    % parameter choice, negative values indicate exceptions.
    if(choice_param_id>0)
      fprintf('Choice was %d', choice_param_id);
      superpixel_param_choices{index}(end+1) = choice_param_id;
      choice_param.method = choose_config.param(choice_param_id).method;
      choice_param.seg_suffix = choose_config.param(choice_param_id).seg_suffix;
      check_for_dir([superpixel_param_choice_dir, image_sub_dir]);
      save2([superpixel_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
        'choice_param');
      
      tile_id = tile_id + 1;
    else
      if(choice_param_id==-1)
        break;
      end
    end
  end
  
  if(choice_param_id~=-1)
    % no exceptions, goto next plane
    index = index + 1;
    tile_id = 1;
    fprintf('\n');
  end
end
fprintf('done.\n');

return
end
