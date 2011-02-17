function exp_segment_cascade_medulla_misc()

c1{1,1} = 'num_b<'; c1{1,2} = '20 1 40';
c1{2,1} = 'num_b>'; c1{2,2} = '100 0 40';

c2{1,1} = 'num_b<'; c2{1,2} = '10 1 5';
c2{2,1} = 'num_b>'; c2{2,2} = '100 0 40';

c3{1,1} = 'num_b<'; c3{1,2} = '30 1 30';
c3{2,1} = 'num_b>'; c3{2,2} = '50 0 20';

c4{1,1} = 'num_b<'; c4{1,2} = '10 1 1';
c4{2,1} = 'num_b>'; c4{2,2} = '50 0 20';

c5{1,1} = 'num_b<'; c5{1,2} = '30 1 100';

cm{1,1} = 'merge_along_longest_boundary'; cm{1,2} = '10 0.5 20';

for i = 12:16
  sp = load(sprintf('../../medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/grayscale_ladder/a.%03d.gs_l_sp_T1e-05_L300.gs_amb_T0.08_0.07_b0.BEL776750629509_40.mat', i));
  sp.label_map = uint32(sp.label_map);
  
  b = load(sprintf('../../medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/BEL/a.%03d.BEL.c_tbb5_mf7.mat', i));
  b.boundary = uint8(255*b.boundary); b.boundary = 255 - b.boundary;
  
  m = load(sprintf('../../medulla.HPF.Leginon.3500x.zhiyuan.fall2008/mitochondria/a.%03d.mitochondria_det_conf.mat', i));
  m.mitochondria_detection_confidence(m.mitochondria_detection_confidence<0) = 0;
  m.mitochondria_detection_confidence = uint8(10*m.mitochondria_detection_confidence);
  image = imread(sprintf('/groups/chklovskii/medulla/medulla.HPF.Leginon.3500x.zhiyuan.fall2008/a.%03d.tif',i));
  
  l1 = remove_merged_boundaries_2D(segment_cascade(sp.label_map, b.boundary, c1));
  l2 = remove_merged_boundaries_2D(segment_cascade(l1, b.boundary, c2));
  l3 = remove_merged_boundaries_2D(segment_cascade(l2, b.boundary, c3));
  l4 = remove_merged_boundaries_2D(segment_cascade(l3, b.boundary, c4));
  l5 = remove_merged_boundaries_2D(segment_cascade(l4, b.boundary, c5));
  
  lf.label_map = remove_merged_boundaries_2D(merge_mitochondria_to_mother_cell(...
    l5, b.boundary, m.mitochondria_detection_confidence, cm));
  
  lf.label_map = double(lf.label_map);
  save(sprintf('../../medulla.HPF.Leginon.3500x.zhiyuan.fall2008/2D_segmentation_results/segment_cascade/a.%03d.seg_cascade_exp.mat', i), '-STRUCT', 'lf');
  
  [py, px] = find(lf.label_map==0);
  figure; imshow(image);
  hold on; plot(px, py, '.'); hold off;
end