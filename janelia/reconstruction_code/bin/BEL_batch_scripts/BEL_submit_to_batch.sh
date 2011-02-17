#!/bin/csh

echo "$EM_CODE_DIR/bin/BEL_batch_scripts/make_sub_dirs.sh $argv[2-]"
$EM_CODE_DIR/bin/BEL_batch_scripts/make_sub_dirs.sh $argv[2-]

$EM_CODE_DIR/bin/BEL_batch_scripts/generate_BEL_batch.sh $argv[1-] > b.sh

chmod u+x b.sh

./b.sh

