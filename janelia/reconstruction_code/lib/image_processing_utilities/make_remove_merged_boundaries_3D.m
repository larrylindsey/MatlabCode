function make_remove_merged_boundaries_3D()
mex -output remove_merged_boundaries_3D ...
    -I./ ...
    -L./ ...
    CFLAGS="\$CFLAGS -O3" ...
    remove_merged_boundaries_3D.mex.cpp
return;
end