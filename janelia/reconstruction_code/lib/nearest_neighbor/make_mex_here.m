function make_mex_here()

mex get_knn_distance_L1_c.cpp
mex get_nn_distance_L1_c.cpp
mex get_nn_distance_L1_weighted_c.cpp
mex -output pdist2 pdist2.mex.cpp

return
end
