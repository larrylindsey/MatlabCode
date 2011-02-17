function make_mex_here()

mex compute_segmentation_from_superpixels_with_min_area_c.cpp

mex -output compute_segmentation_hierarchy_from_watershed_with_min_area_c ...
  -I../lib/merge_sets ...
  CFLAGS="\$CFLAGS -O3" ...
  compute_segmentation_hierarchy_from_watershed_with_min_area_c.cpp

mex seg_frm_superpixel_boundary_vs_interior_FisherLD.cpp
mex segmentation_from_superpixels_mean_boundary_value.cpp
mex segmentation_from_superpixels_median_boundary_value.cpp
mex segment_mean_boundary_value.cpp

make_segmentation_2D_graphcut

make_segmentation_from_superpixels_mean_boundary_value_RAG

mex -output get_rand_score2 ...
  CFLAGS="\$CFLAGS -O3" ...
  get_rand_score2.mex.cpp

mex -output get_normalized_rand_score ...
  CFLAGS="\$CFLAGS -O3" ...
  -I../lib/hash_functions ...
  get_normalized_rand_score.mex.cpp

return
end
