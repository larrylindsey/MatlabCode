  function make_mex_here()

mex -output shrink_segment_levelset_ChanVese ...
  -I../../lib/hash_functions ...
  -I../../lib/levelsets ...
  -L../../lib/levelsets ...
  CFLAGS="\$CFLAGS -O3" ...
  -llevelsets ...
  shrink_segment_levelset_ChanVese.mex.cpp

return
end