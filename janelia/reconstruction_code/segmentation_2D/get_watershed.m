function watershed_img_c = get_watershed(watershed_file_name, boundary_for_watershed)
  if(exist(watershed_file_name, 'file')==2)
    load2(watershed_file_name,  'watershed_img_c');
  else
    watershed_img_c = watershed(boundary_for_watershed, 4);
    save2(watershed_file_name,  'watershed_img_c'); % direct grayscale
  end;
  return;
end

