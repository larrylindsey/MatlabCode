function roi = get_roi_uni_tile(config)
% roi = get_roi_uni_tile(config)
% Compute a "maximalish" square overlap region within slices of an aligned
% stack. There is one tile per section (non trakEM xml)
% roi.xmin, xmax ymin ymax
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  04112008  init code
%

global config_global

if(config.is_verbose)
  fprintf('START: get_roi_uni_tile\n');
end

stack_config = config.stack;

image_dir = get_stack_dir(config);

if(isfield(stack_config.align, 'is_prealigned') && ~stack_config.align.is_prealigned)
  align_config = stack_config.align;
  align_xg_file_name = [image_dir, align_config.xg_file_name];
  xg = read_xf(align_xg_file_name);
  save([image_dir, align_config.xg_save_name], 'xg');
  slice_mask = [];
  for i = 1:length(stack_config.case_ids)
    case_id = stack_config.case_ids(i);
    img = im2double(imread(sprintf([image_dir, stack_config.image_prefix, ...
      stack_config.image_suffix], case_id)));
    margin.left = 1-align_config.margin;
    margin.right = size(img,2)+align_config.margin;
    margin.top = 1-align_config.margin;
    margin.bottom = size(img,1)+align_config.margin;
    s_m = ones(size(img));
    if(isempty(slice_mask))
      align_tform = xf_2_tform(xg(i,:), size(img, 1), size(img, 2));
      slice_mask = apply_tform(s_m, align_tform, margin, 'bilinear')>0;
    else
      align_tform = xf_2_tform(xg(i,:), size(img, 1), size(img, 2));
      slice_mask = slice_mask .* (apply_tform(s_m, align_tform, margin, 'bilinear')>0);
    end;
  end
else
  slice_mask = [];
  for case_id = stack_config.case_ids
    if(config.is_verbose)
      fprintf('%d ', case_id);
      if(mod(case_id, 20)==0)
        fprintf('\n');
      end;
    end

    img = im2double(imread(sprintf([image_dir, stack_config.image_prefix, ...
      stack_config.image_suffix], case_id)));
    if(isempty(slice_mask))
      slice_mask = img>0;
    else
      slice_mask = slice_mask .* (img>0);
    end
  end;
end

if(config.is_verbose_figures)
  figure(1); imshow(slice_mask);
end;

slice_mask_i = 1-slice_mask;
slice_mask_i = bwareaopen(slice_mask_i, 10000);

slice_mask_d = bwdist(slice_mask_i);

[radius, c_index] = max(slice_mask_d(:));

if(radius==inf)
  roi.xmin = 1;
  roi.xmax = size(slice_mask,2);
  roi.ymin = 1;
  roi.ymax = size(slice_mask,1);
  
  if(config.is_verbose)
    fprintf('STOP: get_roi_uni_tile\n');
  end
  return;
end

c_x = floor((c_index-1)/size(slice_mask,1))+1;
c_y = mod((c_index-1), size(slice_mask,1))+1;

width_max = radius-1;
width_min = floor(radius/sqrt(2))-1;

width = floor((width_max+width_min)/2);
while(abs(width_min-width_max)>5)
  roi.xmin = c_x-width;
  roi.xmax = c_x+width;
  roi.ymin = c_y-width;
  roi.ymax = c_y+width;
  roi_mask = zeros(size(slice_mask));
  roi_mask(roi.ymin:roi.ymax, roi.xmin:roi.xmax) = 1;
  if(nnz(roi_mask==1 & slice_mask==0)==0)
    width_min = width;
    width = floor((width+width_max)/2);
  else
    width_max = width;
    width = floor((width_min+width)/2);
  end
end

roi.xmin = c_x-width;
roi.xmax = c_x+width;
roi.ymin = c_y-width;
roi.ymax = c_y+width;

if(config.is_verbose_figures)
  roi_mask = zeros(size(slice_mask));
  roi_mask(roi.ymin:roi.ymax, roi.xmin:roi.xmax) = 1;
  figure(2); imshow(roi_mask);
end;

if(config.is_verbose)
  fprintf('STOP: get_roi_uni_tile\n');
end

return

end
