#!/bin/csh
foreach i ( $argv[1-] )
 echo $i
 cd $i
 ls -l data/*.tif | awk -f $EM_CODE_DIR/bin/BEL_cluster_scripts/w.awk | xargs -i\{\} mkdir data/\{\}
 ls -l data/*.tif | awk -f $EM_CODE_DIR/bin/BEL_cluster_scripts/w.awk | xargs -i\{\} mkdir data/\{\}/canny1
 ls -l data/*.tif | awk -f $EM_CODE_DIR/bin/BEL_cluster_scripts/w.awk | xargs -i\{\} mkdir data/\{\}/gray
 cd ..
 end
