function superpixel_param_choices = gui_choose_segmentation_param(config)
% superpixel_param_choices = gui_choose_segmentation_param(config)
% GUI for choosing segmentation parameters from a user-defined set
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  111908  init code: user specified binary decision tree
%

stack_config = config.stack;
seg_config = config.segmentation_2D;
choose_config = config.segmentation_choose;
tree_config = choose_config.decision_tree;

prev_dir = pwd2;
cd(get_reconstruction_dir(config));
if(exist(choose_config.dir, 'dir')~=7)
  mkdir2(choose_config.dir);
end;
cd(prev_dir);

seg_dir = [get_reconstruction_dir(config), seg_config.dir];

MARGIN = 0.005;
MAX_WIDTH = 1/tree_config.type;

superpixel_param_choices = cell(1,length(stack_config.case_ids));
for index = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(index);
  fprintf('%d:\n', case_id);
  
  [images, image_prefixes] = get_image_from_stack(config, case_id);
  
  for tile_id = 1:length(images)
    fprintf('tile %d ', tile_id);
    
    image = images{tile_id};
    image_prefix = image_prefixes{tile_id};
    
    image = imresize(image, 1/stack_config.segmentation_scale, ...
      stack_config.segmentation_scaling_method);
    display_width = MAX_WIDTH - 2*MARGIN;
    
    current_node_id = 1;
    while(current_node_id>0)
      h_fig = figure(1);
      clf;
      set(h_fig, 'Name', ...
        ['Image, plane: ', num2str(case_id), ' tile ', num2str(tile_id)]);
      for i = 1:tree_config.type
        param_id = tree_config.node(current_node_id).param_id(i);
        seg = load([seg_dir, choose_config.param(param_id).method, '/', ...
          image_prefix, choose_config.param(param_id).seg_suffix, '.mat']);
        subplot('position', [MARGIN + MAX_WIDTH*(i-1), MARGIN, display_width, 1]);
        imshow(image);
        [py, px] = find(seg.label_map==0);
        hold on; plot(px, py, '.'); hold off;
        colormap('gray');
      end
      user_choice = 0;
      while(user_choice==0)
        press_type = waitforbuttonpress;
        while(press_type==0)
          press_type = waitforbuttonpress;
        end
        key_press = get(h_fig, 'CurrentCharacter');
        for k = 1:length(choose_config.input_keys)
          if(strcmp(key_press, choose_config.input_keys{k})==1)
            user_choice = k;
            break;
          end
        end
      end
      current_node_id = tree_config.node(current_node_id).node_id(user_choice);
    end
    choice_param_id = -current_node_id;
    
    fprintf('Choice was %d', choice_param_id);
    superpixel_param_choices{index}(end+1) = choice_param_id;
  end
  fprintf('\n');
end
fprintf('done.\n');

return
end
