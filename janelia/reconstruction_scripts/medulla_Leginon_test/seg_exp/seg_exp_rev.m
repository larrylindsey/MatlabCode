
load config.medulla.mat

for i = 12:16
  
  image_prefix = sprintf('a.%03d', i);
  fprintf('image_prefix: %s\n', image_prefix);
  
  im = imread(['/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/', image_prefix, '.tif']);
  load(['/groups/chklovskii/medulla/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/mitochondria/', image_prefix, '.mitochondria_det_conf.mat']);
  sg2 = load2(['../../../medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/prune_classify/', image_prefix, ...
    '.apc_sp_region_max_3_0.35', ...
    '.gs_l_sp_T0.25_L300', ...
    '.gs_amdb_sp_T0.35_0.34_b5', ...
    '.gs_amb_T0.25_0.24_b0', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg', ...
    '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg.mat']);
  
  % f = filter_image2(im, 'db_nwo700_cbf2_LDA7_mf7_ps1\_0.25_abf_neg', config);
  % save([image_prefix, '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg.mat'], 'f');
  load([image_prefix, '.db_nwo700_cbf2_LDA7_mf7_ps1_0.25_abf_neg.mat']);
  
  b = 1-f;
  
  m = mitochondria_detection_confidence>0;
  m = imerode(m, strel('disk', 5));
  m = bwareaopen(m, 40);
  
  s = b<0.3;
  s(m>0) = 0;
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map_0 = double(remove_merged_boundaries_2D(uint32(label_map)));
  fprintf('number of segments: %d\n', length(unique(label_map_0(:)))-1);

  s = b<0.3;
  s(m>0) = 0;
  s(1:end-1, 1:end-1) = s(2:end, 2:end);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n2 = segment_and(label_map_0, label_map, b<0.5, 150);
  
  s = b<0.3;
  s(m>0) = 0;
  s(2:end, 2:end) = s(1:end-1, 1:end-1);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n3 = segment_and(label_map_n2, label_map, b<0.5, 150);

  s = b<0.3;
  s(m>0) = 0;
  s(2:end, 1:end-1) = s(1:end-1, 2:end);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n4 = segment_and(label_map_n3, label_map, b<0.5, 150);

  s = b<0.3;
  s(m>0) = 0;
  s(1:end-1, 2:end) = s(2:end, 1:end-1);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n5 = segment_and(label_map_n4, label_map, b<0.5, 150);

  s = b<0.3;
  s(m>0) = 0;
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 150);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n6 = segment_and(label_map_n5, label_map, b<0.5, 150);

  s = b<0.3;
  s(m>0) = 0;
  s(1:end-1, 2:end) = s(1:end-1, 1:end-1);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n7 = segment_and(label_map_n6, label_map, b<0.5, 150);

  s = b<0.3;
  s(m>0) = 0;
  s(1:end-1, 2:end) = s(2:end, 2:end);
  
  sde = imerode(imdilate(s, strel('disk', 2)), strel('disk', 2));
  
  sdea = bwareaopen(sde, 20);
  
  bs = imimposemin(b, sdea, 4);
  
  ws_bs = watershed(bs, 4);
  label_map = compute_segmentation_hierarchy_from_watershed_with_min_area_c(...
    ws_bs, bs, 0.1, 200);
  label_map = double(remove_merged_boundaries_2D(uint32(label_map)));

  label_map_n8 = segment_and(label_map_n7, label_map, b<0.5, 150);

  label_map_n9 = segment_and(label_map_n8, sg2.label_map, b<0.5, 150);
  
%  plot_segment_boundaries(im, label_map_n7);
  
  label_map = label_map_n9;
  fprintf('number of segments after mergers: %d\n', length(unique(label_map(:)))-1);
  save(['~/research/em_reconstruction_pipeline/reconstructions/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/seeded_watershed/', image_prefix, ...
    '.sws_0.3_2_20_4_200_and_APC2_0.35_nm_rev.mat'], 'label_map');
end