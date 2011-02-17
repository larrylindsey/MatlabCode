function adjust_imported_proofread_superpixel_map(config)

if(isfield(config.stack, 'image_structure') && ~isempty(config.stack.image_structure))
  adjust_imported_proofread_superpixel_map_tile_patches(config);
else
  import_proofread_superpixel_map_uni_tile(config);
end

return
end
