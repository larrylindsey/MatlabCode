function generate_al_uni_tile_coarse(config)
% generate_al_uni_tile_coarse(config)
% generate coarse al when there is one tile per image (non trakEM xml)
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

fprintf('Constructing coarse al without roi for preview ...\n');
fprintf('NOTE: run 810.35 to produce final al for proofreader.\n');

stack_config = config.stack;

image_dir = get_stack_dir(config);

% stack is not prealigned then compute an ROI for the aligned stack if
% needed
if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned ...
    && isempty(stack_config.align.roi))
  align_roi_file_name = [get_region_dir(config), config.stack.align.roi_file_name];
  load2(align_roi_file_name, 'align_roi');
end

case_id = stack_config.case_ids(1);
fprintf('%d ', case_id);

image_1 = get_image_from_stack(config, case_id);
image_1 = uint8(image_1{1}*255); % just take the first tile. Temporary fix - must be corrected to load multiple tiles

if(~isempty(stack_config.roi))
  image_1 = image_1(stack_config.roi.ymin:stack_config.roi.ymax, ...
    stack_config.roi.xmin:stack_config.roi.xmax);
end;

if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xg_file_name = [image_dir, align_config.xg_file_name];
  xg = read_xf(align_xg_file_name);
  save2([image_dir, align_config.xg_save_name], 'xg');
  
  margin.left = 1-align_config.margin;
  margin.right = size(image_1,2)+align_config.margin;
  margin.top = 1-align_config.margin;
  margin.bottom = size(image_1,1)+align_config.margin;

  align_tform = xf_2_tform(xg(1,:), size(image_1, 1), size(image_1, 2));
  image_1 = uint8(apply_tform(image_1, align_tform, margin, 'bilinear'));
  image_1 = image_1(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
end;

save_file_name = [get_region_dir(config), config.coarse_al.file_name];
if(exist(save_file_name, 'file')==2)
  delete(save_file_name);
end
al = image_1;
al = imresize(al, 1/config.coarse_al.scale, 'bilinear');
imwrite(al, save_file_name, 'Compression', 'lzw', 'WriteMode', 'append');

for i = 2:length(stack_config.case_ids)
  case_id = stack_config.case_ids(i);
  fprintf('%d ', case_id);
  if(mod(i,20)==0)
    fprintf('\n');
  end;
  
  image_2 = get_image_from_stack(config, case_id);
  image_2 = uint8(image_2{1}*255); % just take the first tile. Temporary fix - must be corrected to load multiple tiles
  
  if(~isempty(stack_config.roi))
    image_2 = image_2(stack_config.roi.ymin:stack_config.roi.ymax, ...
      stack_config.roi.xmin:stack_config.roi.xmax);
  end;
  
  if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
    margin.left = 1-align_config.margin;
    margin.right = size(image_2,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(image_2,1)+align_config.margin;

    align_tform = xf_2_tform(xg(i,:), size(image_2, 1), size(image_2, 2));
    image_2 = uint8(apply_tform(image_2, align_tform, margin, 'bilinear'));
    image_2 = image_2(align_roi.ymin:align_roi.ymax, align_roi.xmin:align_roi.xmax);
  end;
  
  al = image_2;
  al = imresize(al, 1/config.coarse_al.scale, 'bilinear');
  imwrite(al, save_file_name, 'Compression', 'lzw', 'WriteMode', 'append');
end;

output_to_proofreader_stitched_grayscale(al, config);

fprintf('\ndone.\n');

