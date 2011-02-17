function output_3D_segmentation_to_raveler( ambc_thresh )

fprintf('START: output_3D_segmentation_to_raveler\n');

scale = 1;

% output_dir = '~/temp/FIB.46.162.amb.T120/';
% stack_dir = '/groups/chklovskii/chklovskiilab/electron_microscopy_data/rat_cortex.FIB.013009/';
% image_prefix = 'FIBSLICE%04d';
% image_suffix = '.TIF';
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
% % seg_seg_dir = [seg_dir, 'stitch/'];
% % seg_prefix = 'stitch_seg.%d.%d.mf5.ws.T60.ld.T60_L1000.46_162';

% output_dir = '~/temp/fly_larva.FIB.D01_5x5x20nm_sn0309_73.061809.raveler/';
% stack_dir = '/groups/chklovskii/chklovskiilab/electron_microscopy_data/fly_larva.FIB.D01_5x5x20nm_sn0309_73.061809/';
% image_prefix = 'a.%03d';
% roi.xmin = 201; roi.xmax = 1800; roi.ymin = 201; roi.ymax = 1800;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = [1 30 59]; end_indexs = [30 59 84];
% % start_indexs = 1; end_indexs = 30;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/fly_larva.FIB.D01_5x5x20nm_sn0309_73.061809/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf5_heq.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% % seg_prefix = 'seg_stack.%d.%d.mf5_heq.ws.T1.ld.T1_L1000.amb.T140.ambcs.T200.1_84';
% seg_prefix = 'seg_stack.%d.%d.mf5_heq.ws.T1.ambl.T160.L50.A1000.ambcs.T200.1_84';

% % % output_dir = '/media/FreeAgent_Drive/temp/fly_larva.FIB.5x5x50nm/';
% % output_dir = '~/temp/fly_larva.FIB.5x5x50nm/';
% % % stack_dir = '/media/FreeAgent_Drive/em_reconstruction/data_em_images/fly_larva.FIB.5x5x5nm/';
% % stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/fly_larva.FIB.5x5x5nm/';
% % image_prefix = 'dn8.%03d';
% % % roi.xmin = 201; roi.xmax = 1800; roi.ymin = 201; roi.ymax = 1800;
% % filter_version = 'heq';
% % image_suffix = '.tif';
% % start_indexs = [1 30 59 88 117 146 175];
% % end_indexs = [30 59 88 117 146 175 204];
% % % start_indexs = [1 30 59];
% % % end_indexs = [30 59 88];
% % % seg_dir = '/media/FreeAgent_Drive/em_reconstruction/reconstructions/fly_larva.FIB.5x5x5nm/3D_segmentation_results/';
% % seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/fly_larva.FIB.5x5x5nm/3D_segmentation_results/';
% % seg_ws_dir = [seg_dir, 'watershed/'];
% % seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.neg.ws.T1.raw'];
% % seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% % % seg_prefix = 'seg_stack.%d.%d.neg.ws.T1.ambl.T160.L50.A1000';
% % % seg_prefix = 'seg_stack.%d.%d.neg.ws.T1.ambl.T160.L50.A1000.ambcs.T200.1_88';
% % seg_prefix = 'seg_stack.%d.%d.neg.ws.T1.ambl.T160.L50.A1000.ambc.T160.1_204';

% output_dir = '~/temp/tomo_pilot_11112009a/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_pilot_11112009a/';
% image_prefix = 'a.%03d';
% roi.xmin = 100; roi.xmax = 350; roi.ymin = 150; roi.ymax = 400;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% % start_indexs = [1 9];
% % end_indexs = [8 17];
% start_indexs = [1 9 18 27 36 45];
% end_indexs = [8 17 26 35 44 53];
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_pilot_11112009a/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf3.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambc.T98.1_53';

% output_dir = '~/temp/tomo_FIB_recon_ambcs/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_FIB_recon/';
% image_prefix = 'r.%03d';
% roi.xmin = 25; roi.xmax = 225; roi.ymin = 25; roi.ymax = 225;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = 10;
% end_indexs = 261;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_FIB_recon/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf3.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% % seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambc.T140.10_261';
% seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambcs.T140.10_261';

% output_dir = '~/temp/tomo_fly_larva_FIB_11172009_ambc/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_fly_larva_FIB_11172009/';
% image_prefix = 'r.%03d';
% roi.xmin = 25; roi.xmax = 275; roi.ymin = 25; roi.ymax = 275;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = 10;
% end_indexs = 306;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_fly_larva_FIB_11172009/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf3.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambc.T130.10_306';
% % seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambcs.T140.10_306';

% output_dir = '~/temp/tomo_fly_larva_FIB_orig_11172009_ambc/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_fly_larva_FIB_orig_11172009/';
% image_prefix = 'o.%03d';
% roi.xmin = 25; roi.xmax = 275; roi.ymin = 25; roi.ymax = 275;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = 10;
% end_indexs = 306;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_fly_larva_FIB_orig_11172009/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.mf3.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambc.T130.10_306';
% % seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambcs.T140.10_306';

% output_dir = '~/temp/tomo_fly_larva_30sec_smoothed_11182009_ambc.190/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_fly_larva_30sec_smoothed_11182009/';
% image_prefix = 's.%03d';
% roi.xmin = 25; roi.xmax = 350; roi.ymin = 25; roi.ymax = 350;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = [1 10 19 28 37 46 55 64 73 82 91 100 109 118 127 136 145 154 163 172 181 190 199 208 217 226 235 244 253 262];
% end_indexs = [9 18 27 36 45 54 63 72 81 90 99 108 117 126 135 144 153 162 171 180 189 198 207 216 225 234 243 252 261 270];
% % start_indexs = [1 10 19 28 37 46 55 64 73 82];
% % end_indexs = [9 18 27 36 45 54 63 72 81 90];
% % start_indexs = [1 10];
% % end_indexs = [9 18];
% % start_indexs = 10;
% % end_indexs = 18;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_fly_larva_30sec_smoothed_11182009/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.heq_mf3.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.heq_mf3.ws.T1.ambl.T170.L50.A1000.ambc.T190.1_270';
% % seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambcs.T140.10_306';

% output_dir = '~/temp/tomo_fly_larva_orig_images_1/';
% stack_dir = '~/research/em_reconstruction_pipeline/data_em_images/tomo_fly_larva_orig_images_1/';
% image_prefix = 'o.%03d';
% roi.xmin = 25; roi.xmax = 350; roi.ymin = 25; roi.ymax = 350;
% filter_version = 'heq_neg';
% image_suffix = '.tif';
% start_indexs = 1;
% end_indexs = 30;
% seg_dir = '/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/reconstructions/tomo_fly_larva_orig_images_1/3D_segmentation_results/';
% seg_ws_dir = [seg_dir, 'watershed/'];
% seg_ws_name = [seg_ws_dir, 'seg_stack.%d.%d.heq_d4_mf8.ws.T1.raw'];
% seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
% seg_prefix = 'seg_stack.%d.%d.heq_d4_mf8.ws.T1.ambl.T120.L50.A1000.ambc.T200.1_30';
% % seg_prefix = 'seg_stack.%d.%d.mf3.ws.T1.ambl.T120.L50.A1000.ambcs.T140.10_306';

if nargin < 1,
    ambc_thresh = 130;
end

output_dir = ['~/em_reconstructions/reconstructions/spams_10/for_proofreading/T' sprintf('%d',ambc_thresh) 'adj/'];
system(['mkdir -p ' output_dir]);
stack_dir = '/groups/chklovskii/home/nuneziglesiasj/Projects/em_denoising/data/10x10x10_cropped/';
image_prefix = 'a.%03d';
roi = [];
filter_version = '';
image_suffix = '.tif';
start_indexs = 0;
end_indexs = 499;
stack_filter_version = 'dict10x10x10.iter1085.mf3';
seg_dir = '/groups/chklovskii/home/nuneziglesiasj/em_reconstructions/reconstructions/spams_10/3D_segmentation_results/';
seg_ws_dir = [seg_dir, 'watershed/'];
seg_ws_name = [seg_ws_dir, 'seg_stack.0.499.', stack_filter_version, '.ws.T10.raw'];
seg_seg_dir = [seg_dir, 'agglo_mean_boundary/'];
sp_prefix = ['seg_stack.0.499.', stack_filter_version, '.ws.T10.ambl.T90.L50.A1000'];
seg_prefix = ['seg_stack.0.499.', stack_filter_version, '.ws.T10.ambl.T90.L50.A1000.ambc.T' sprintf('%d',ambc_thresh) '.0_499'];

sp_seg_name = [seg_seg_dir, sp_prefix, '.raw'];
seg_seg_name = [seg_seg_dir, seg_prefix, '.raw'];

check_for_dir(output_dir);
system(['rm -rf ', output_dir, '*']);

fprintf('Dumping grayscale maps ...\n');
grayscale_dir = [output_dir, 'grayscale_maps/'];
mkdir(grayscale_dir);
for i = start_indexs(1):end_indexs(end)
  fprintf('i: %d\n', i);
  image = im2double(imread(sprintf([stack_dir, image_prefix, image_suffix], i)));
  if(~isempty(roi))
    image = im2double(image(roi.ymin:roi.ymax, roi.xmin:roi.xmax));
  end
  image = image / max(image(:));
  image = filter_image2(image, filter_version);
  if(scale~=1)
    image = imresize(image, 1/scale, 'bilinear');
  end
  imwrite(image, sprintf([grayscale_dir, 'image.v1.%05d.png'], i));
end

% segment maps
fprintf('Dumping superpixel and segment maps ...\n');
superpixel_2_segment_map = [];
segment_2_body_map = [0 0];
superpixel_dir = [output_dir, 'superpixel_maps/'];
mkdir(superpixel_dir);
segment_offset = 0;
for ss_i = 1:length(start_indexs)
  start_index = start_indexs(ss_i);
  end_index = end_indexs(ss_i);
  [seg_ws, err] = fread_raw_array_mex(...
    sprintf(seg_ws_name, start_index, end_index));
  if(err==1)
    fprintf('Could not read watershed file:%s\n', seg_ws_name);
    error('foo');
  end

  if(~isempty(sp_seg_name))
    [ws_2_sp, err] = fread_raw_array_mex(...
      sprintf(sp_seg_name, start_index, end_index));
    if(err==1)
      fprintf('Could not read superpixel file:%s\n', sp_seg_name);
      error('foo');
    end
  else
    ws_2_sp = [];
  end
  
  [ws_2_seg, err] = fread_raw_array_mex(...
    sprintf(seg_seg_name, start_index, end_index));
  if(err==1)
    fprintf('Could not read segment file:%s\n', seg_seg_name);
    error('foo');  
  end
  
  if(ss_i~=1)
    start_index_0 = start_index;% + 1;
  else
    start_index_0 = start_index;
  end
  for i = start_index_0:end_index
    fprintf('i: %d\n', i);
    ws_z_3D = seg_ws(:,:, i-start_index+1)';

    if(scale~=1)
      ws_z_3D = imresize(ws_z_3D, 1/scale, 'nearest');
    end
    if(isempty(ws_2_sp))
      sp_z_2D = relabel_connected_components_2D(uint32(ws_z_3D));
    else
      sp_z_2D = relabel_connected_components_2D(uint32(ws_2_sp(1+ws_z_3D)));
    end
    
    sp_prop = regionprops(sp_z_2D);
    to_delete_sp_ids = find([sp_prop(:).Area]<20);
    sp_z_2D(ismember(sp_z_2D, to_delete_sp_ids)) = 0;
    
    body_z_3D = ws_2_seg(1+ws_z_3D);
    seg_z_2D = relabel_connected_components_2D(uint32(body_z_3D));
    
    seg_z_2D(seg_z_2D>0) = seg_z_2D(seg_z_2D>0) + segment_offset;
    segment_offset = max(segment_offset, max(seg_z_2D(:)));
    
    sp_b = draw_segment_boundaries_c(double(sp_z_2D));
    sp_z_2D(sp_b==1) = 0;

    sp_2_seg_z_2D = unique([sp_z_2D(:), seg_z_2D(:)], 'rows');
    sp_2_seg_z_2D = sp_2_seg_z_2D(min(sp_2_seg_z_2D, [], 2)>0, :);
    sp_2_seg_z_2D = [0 0; sp_2_seg_z_2D]; %#ok<AGROW>
    
    seg_2_body_z_2D = unique([seg_z_2D(:), body_z_3D(:)], 'rows');
    seg_2_body_z_2D = seg_2_body_z_2D(min(seg_2_body_z_2D, [], 2)>0, :);
    seg_2_body_z_2D = [0 0; seg_2_body_z_2D]; %#ok<AGROW>

    [~,f] = unique(sp_2_seg_z_2D(:,1), 'first');
    [~,l] = unique(sp_2_seg_z_2D(:,1), 'last');
    nnz(f~=l)
    [~,f] = unique(seg_2_body_z_2D(:,1), 'first');
    [~,l] = unique(seg_2_body_z_2D(:,1), 'last');
    nnz(f~=l)
    
    imwrite(uint16(sp_z_2D), ...
      sprintf([superpixel_dir, 'superpixel_map.v1.%05d.png'], i), ...
      'BitDepth', 16);
    
    superpixel_2_segment_map = [superpixel_2_segment_map;
      repmat(i, [size(sp_2_seg_z_2D,1), 1]), sp_2_seg_z_2D]; %#ok<AGROW>    
    
    segment_2_body_map = [segment_2_body_map; seg_2_body_z_2D]; %#ok<AGROW>
  end
end

segment_2_body_map = unique(segment_2_body_map, 'rows');

fout = fopen([output_dir, 'superpixel_to_segment_map.txt'], 'wt');
fprintf(fout, '%d %d %d\n', superpixel_2_segment_map');
fclose(fout);

fout = fopen([output_dir, 'segment_to_body_map.txt'], 'wt');
fprintf(fout, '%d %d\n', segment_2_body_map');
fclose(fout);

fprintf('STOP: output_3D_segmentation_to_raveler\n');
return
end
