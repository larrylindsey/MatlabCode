function display_aligned_SIFT_affine_tile_pair_inter_plane(config)
% display_aligned_SIFT_affine_tile_pair_inter_plane(config)
% Display aligned tile images using transformations computed for each
% plane/section separately. This is part of two stage alignment
%
% Shiv N. Vitaladevuni
% Janelia Farm Research Campus, HHMI
%
% v0  08312008  init. code
%

stack_config = config.stack;
sift_config = config.align.linkage_align.SIFT;
sift_dir = get_sift_dir(config);

for layer_id = 1:length(stack_config.case_ids)-1
  case_id_1 = stack_config.case_ids(layer_id);
  case_id_2 = stack_config.case_ids(layer_id+1);

  [images_1, image_prefixes_1] = get_image_from_stack(config, case_id_1);
  [images_2, image_prefixes_2] = get_image_from_stack(config, case_id_2);

  for tile_1 = 1:length(images_1)
    if(isempty(images_1{tile_1}))
      continue;
    end
    if(sift_config.is_inverted_display)
      images_1{tile_1} = 1-images_1{tile_1};
    end
    for tile_2 = 1:length(images_2)
      if(isempty(images_2{tile_2}))
        continue;
      end
      file_name_suffix = get_file_name_from_tuple(sift_dir, ...
        image_prefixes_1{tile_1}, image_prefixes_2{tile_2}, 'st.');
      file_name = [file_name_suffix, '.mat'];
%       file_name=[sift_dir, 'sift_transform.tilepair.', num2str(case_id_1), ...
%         '_', num2str(case_id_2),'.mat'];
      try
        load2(file_name, 'transforms_tp');
      catch
        continue;
      end
      if(sift_config.is_inverted_display)
        image_2 = 1-images_2{tile_2};
      else
        image_2 = images_2{tile_2};
      end
      xdata = {};
      ydata = {};
      tform_1 = maketform('affine', ...
        reshape(transforms_tp(1,:), [2 3])');
      [image_t_1, xdata{1}, ydata{1}] = imtransform(images_1{tile_1}, tform_1);
      minx = min(xdata{1}); miny = min(ydata{1});
      maxx = max(xdata{1}); maxy = max(ydata{1});
      
      tform_2 = maketform('affine', ...
        reshape(transforms_tp(2,:), [2 3])');
      [image_t_2, xdata{2}, ydata{2}] = imtransform(image_2, tform_2);
      minx = min([minx, xdata{2}]);
      miny = min([miny, ydata{2}]);
      maxx = max([maxx, xdata{2}]);
      maxy = max([maxy, ydata{2}]);
      for t = 1:2
        xdata{t} = round(xdata{t}-minx+1);
        ydata{t} = round(ydata{t}-miny+1);
      end
      maxx = round(maxx - minx + 1);
      maxy = round(maxy - miny + 1);
      
      tile_image = zeros(maxy, maxx, 3);
      tile_image(ydata{1}(1):ydata{1}(2), xdata{1}(1):xdata{1}(2), 1) = ...
        image_t_1;
      tile_image(ydata{2}(1):ydata{2}(2), xdata{2}(1):xdata{2}(2), 2) = ...
        image_t_2;
      figure; imshow(tile_image);
    end
  end
end

return
end
