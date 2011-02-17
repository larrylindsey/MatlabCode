function generate_matlab_global_tforms_multi_tile_patches_compile(config)
% generate_matlab_global_tforms_multi_tile_patches(config)
% generate MATLAB style transforms for all the patches when there are
% multiple tiles and patches per section (trakEM xml). These transforms
% will be used for rendering al and cat.
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
% v1  09302008  modified to multiple tile per section (trakEM xml)
%

if(config.is_verbose)
  fprintf('START: generate_matlab_global_tforms_multi_tile_patches_compile\n');
end

stack_config = config.stack;

minx = inf;
miny = inf;
maxx = -inf;
maxy = -inf;
tforms_all = {};
xdata_all = {};
ydata_all = {};
for layer_id = 1:length(stack_config.case_ids)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  
  tx = load2([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.', num2str(case_id), ...
    '.mat'], 'xdata', 'ydata', 'tforms');

  for tile_id = 1:size(tx.xdata,1)
    for patch_id = 1:size(tx.xdata,2)
    if(isempty(tx.xdata{tile_id, patch_id}))
      continue;
    end
    
    tforms_all{layer_id, tile_id, patch_id} = tx.tforms{tile_id, patch_id}; %#ok<AGROW>
    xdata_all{layer_id, tile_id, patch_id} = tx.xdata{tile_id, patch_id}; %#ok<AGROW>
    ydata_all{layer_id, tile_id, patch_id} = tx.ydata{tile_id, patch_id}; %#ok<AGROW>
    
    minx = min([minx, xdata_all{layer_id, tile_id, patch_id}]);
    miny = min([miny, ydata_all{layer_id, tile_id, patch_id}]);
    maxx = max([maxx, xdata_all{layer_id, tile_id, patch_id}]);
    maxy = max([maxy, ydata_all{layer_id, tile_id, patch_id}]);
    end
  end
end

if(config.is_verbose)
  fprintf('Saving shifted transformations.\n');
end
copy_save_dir = [get_to_be_proofread_dir(config), 'global_stitching_parameters/'];
check_for_dir(copy_save_dir);
for layer_id = 1:size(xdata_all,1)
  case_id = stack_config.case_ids(layer_id);
  fprintf('plane: %d\n', case_id);
  
  xdata = {};
  ydata = {};
  tforms = {};
  for tile_id = 1:size(xdata_all,2)
    for patch_id = 1:size(xdata_all,3)
      if(isempty(xdata_all{layer_id, tile_id, patch_id}))
        continue;
      end
      xdata{tile_id, patch_id} = ...
        round(xdata_all{layer_id, tile_id, patch_id}-minx+1); %#ok<AGROW>
      ydata{tile_id, patch_id} = ...
        round(ydata_all{layer_id, tile_id, patch_id}-miny+1); %#ok<AGROW>
      tforms{tile_id, patch_id} = tforms_all{layer_id, tile_id, patch_id}; %#ok<AGROW>
    end
  end
  
  save([get_global_stitching_param_dir(config), ...
    'global_stitching_parameters.c.', num2str(case_id), ...
    '.mat'], 'xdata', 'ydata', 'tforms');
  save([copy_save_dir, 'global_stitching_parameters.c.', num2str(case_id), ...
    '.mat'], 'xdata', 'ydata', 'tforms');
end
canvas_width = round(maxx - minx + 1); %#ok<NASGU>
canvas_height = round(maxy - miny + 1); %#ok<NASGU>

save2([get_global_stitching_param_dir(config), 'canvas.mat'], ...
  'canvas_height', 'canvas_width');
save2([copy_save_dir, 'canvas.mat'], 'canvas_height', 'canvas_width');

if(config.is_verbose)
  fprintf('STOP: generate_matlab_global_tforms_multi_tile_patches_compile\n');
end
return;
end

