#!/bin/csh

echo "$EM_CODE_DIR/bin/BEL_cluster_scripts/make_sub_dirs.sh $argv[2-]"
$EM_CODE_DIR/bin/BEL_cluster_scripts/make_sub_dirs.sh $argv[2-]

$EM_CODE_DIR/bin/BEL_cluster_scripts/generate_BEL_cluster.sh $argv[1-] > b.sh

chmod u+x b.sh

./b.sh

