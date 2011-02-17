function save_2D_segmentation_maps_from_3D_segmentation()

% start_indexs = 46:29:161;
% end_indexs = 75:29:162;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/rat_cortex.FIB.013009/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf5.ws.T60.raw'];
% seg_sp_dir = [seg_dir, 'ladder/'];
% seg_sp_name = [seg_sp_dir, 'seg_stack.%d.%d.mf5.ws.T60.ld.T60_L1000.raw'];
% % seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% % seg_prefix = 'seg_stack.%d.%d.mf5.ws.T60.ld.T60_L1000.amb.T140';
% % seg_seg_name = [seg_seg_dir, seg_prefix, '.raw'];
% seg_seg_dir = [seg_dir, 'stitch/'];
% seg_prefix = 'stitch_seg.%d.%d.mf5.ws.T60.ld.T60_L1000.amb.T120.46_162';
% seg_seg_name = [seg_seg_dir, seg_prefix, '.raw'];
% % seg_seg_dir = [seg_dir, 'stitch/'];
% % seg_prefix = 'stitch_seg.%d.%d.mf5.ws.T60.ld.T60_L1000.46_162';
% % seg_seg_name = [seg_seg_dir, seg_prefix, '.raw'];

start_indexs = [1 30 59]; end_indexs = [30 59 84];
% start_indexs = 1; end_indexs = 30;
seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/fly_larva.FIB.D01_5x5x20nm_sn0309_73.061809/3D_segmentation_results/';
seg_ws_dir = [seg_dir, 'watershed/'];
seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf5_heq.ws.T1.raw'];
seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.mf5_heq.ws.T1.ld.T1_L1000.amb.T140.ambcs.T200.1_84';
seg_prefix = 'seg_stack.%d.%d.mf5_heq.ws.T1.ambl.T160.L50.A1000.ambcs.T200.1_84';
seg_seg_name = [seg_seg_dir, seg_prefix, '.raw'];

output_dir = sprintf([seg_dir, '2D_maps/', seg_prefix, '/'], ...
  start_indexs(1), end_indexs(end));
mkdir2(output_dir);

system(['rm -rf ', output_dir, '*']);

% segment maps
for ss_i = 1:length(start_indexs)
  start_index = start_indexs(ss_i);
  end_index = end_indexs(ss_i);
  fprintf(sprintf(['watershed file: ', seg_ws_name, '\n'], ...
    start_index, end_index));
  [seg_ws, err] = fread_raw_array_mex(...
    sprintf(seg_ws_name, start_index, end_index));
  if(err==1)
    error('Could not read watershed file.');
  end
  
  fprintf(sprintf(['segment file: ', seg_seg_name, '\n'], ...
    start_index, end_index));
  [sp_2_seg, err] = fread_raw_array_mex(...
    sprintf(seg_seg_name, start_index, end_index));
  if(err==1)
    error('Could not read segment file.');
  end
  
  seg_sp = sp_2_seg(1+seg_ws);
  
  if(ss_i~=1)
    start_index_0 = start_index + 1;
  else
    start_index_0 = start_index;
  end
  for i = start_index_0:end_index
    fprintf('i: %d\n', i);
    sp = seg_sp(:,:, i-start_index+1)';
    sp_id = nonzeros(unique(sp(:)));
    sp_2_new_sp_map = [];
    sp_2_new_sp_map(1+sp_id) = 1:length(sp_id); %#ok<AGROW>
    sp_new = sp_2_new_sp_map(1+sp);
    imwrite(uint16(sp_new), ...
      sprintf([output_dir, 'superpixel_map.v1.%05d.png'], i), ...
      'BitDepth', 16);
  end
end

return
end
