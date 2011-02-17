#!/bin/csh

set a=$argv[1]
while($a <= $argv[2])
 echo qsub -N c$a -j y -o log -b y -cwd -v LD_LIBRARY_PATH '"./copy_global_patchwise_transforms_from_simple_txt_to_mat_files '$a'"'
 @ a = $a + 1
end
