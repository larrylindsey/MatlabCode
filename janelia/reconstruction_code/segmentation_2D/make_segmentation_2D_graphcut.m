% make graphcut-matlab interface

if(isunix)
  mex -o segmentation_2D_graphcut_2 ...
    -L../lib/MRF2.1   -I../lib/MRF2.1 -lMRF segmentation_2D_graphcut_2.mex.cpp
end

if(ispc)
  mex -output segmentation_2D_graphcut_2 ...
    -L../lib/MRF2.1   -I../lib/MRF2.1 segmentation_2D_graphcut_2.mex.cpp
end
