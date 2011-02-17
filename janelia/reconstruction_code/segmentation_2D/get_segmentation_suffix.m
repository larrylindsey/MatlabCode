function [seg_method, seg_suffix] = ...
  get_segmentation_suffix(config, image_prefix)
% get the superpixel suffix for an image, either specified in the config or
% chosen for each image using GUI and stored.

seg_method = config.segmentation_choose.choice.method;
seg_suffix = config.segmentation_choose.choice.seg_suffix;

if(strcmp(seg_method, '#(SEG_CHOICE)')==1)
  choose_config = config.segmentation_choose;
  segmentation_param_choice_dir = get_segmentation_param_choice_dir(config);
  load2([segmentation_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
    'choice_param');
  seg_method = choice_param.method;
end

if(~isempty(strfind(seg_suffix, '#(SEG_CHOICE)')))
  if(exist('choice_param', 'var')==0)
    choose_config = config.segmentation_choose;
    segmentation_param_choice_dir = get_segmentation_param_choice_dir(config);
    load2([segmentation_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
      'choice_param');
  end
  seg_suffix = strrep(seg_suffix, '#(SEG_CHOICE)', ...
    choice_param.seg_suffix);
end

if(~isempty(strfind(seg_suffix, '#(SP_CHOICE)')))
  choose_config = config.superpixel_choose;
  segmentation_param_choice_dir = get_superpixel_param_choice_dir(config);
  load2([segmentation_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
    'choice_param');
  seg_suffix = strrep(seg_suffix, '#(SP_CHOICE)', choice_param.seg_suffix);
end

return
end
