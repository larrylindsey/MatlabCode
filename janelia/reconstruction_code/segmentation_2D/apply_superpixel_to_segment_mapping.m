function seg_label_map = apply_superpixel_to_segment_mapping(...
  sp_label_map, sp_to_seg_label)

sp_to_seg_map = [];
sp_to_seg_map(1+sp_to_seg_label(:,1)) = sp_to_seg_label(:,2);
sp_label_map(sp_label_map<0) = 0;
seg_label_map = sp_to_seg_map(1+sp_label_map);

return
end
