function [superpixel_method, superpixel_suffixes] = ...
  get_superpixel_suffixes(config, image_prefix)
% get the superpixel suffix for an image, either specified in the config or
% chosen for each image using GUI and stored.

seg_config = config.segmentation_2D;

if(~isempty(seg_config.superpixel_suffix) && ~iscell(seg_config.superpixel_suffix))
  seg_config.superpixel_suffix = {seg_config.superpixel_suffix};
end

superpixel_method = [];
superpixel_suffixes = {};
for i = 1:length(seg_config.superpixel_suffix)
  superpixel_suffixes{i} = seg_config.superpixel_suffix{i};
end

if(isempty(seg_config.superpixel_method))
  superpixel_method = config.superpixel_choose.choice.method;
  superpixel_suffixes{1} = config.superpixel_choose.choice.seg_suffix;
end

if(strcmp(superpixel_method, '#(SP_CHOICE)')==1)
  choose_config = config.superpixel_choose;
  superpixel_param_choice_dir = get_superpixel_param_choice_dir(config);
  load2([superpixel_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
    'choice_param');
  superpixel_method = choice_param.method;
  
elseif(isempty(superpixel_method))
  superpixel_method = seg_config.superpixel_method;
  
end

for i = 1:length(superpixel_suffixes)
  if(~isempty(strfind(superpixel_suffixes{i}, '#(SP_CHOICE)')))
    if(exist('choice_param', 'var')==0)
      choose_config = config.superpixel_choose;
      superpixel_param_choice_dir = get_superpixel_param_choice_dir(config);
      load2([superpixel_param_choice_dir, image_prefix, choose_config.save_suffix, '.mat'], ...
        'choice_param');
    end
    superpixel_suffixes{i} = strrep(superpixel_suffixes{i}, '#(SP_CHOICE)', ...
      choice_param.seg_suffix);
  end
end

return
end
