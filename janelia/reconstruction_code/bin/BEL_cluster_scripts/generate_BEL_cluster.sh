#!/bin/csh
set BEL_CONFIG_FILE_SUFFIX=$1

foreach i ( $argv[2-])
 echo cd $i
 echo qsub -N $i -j y -o /dev/null -b y -cwd '"/usr/local/wine/bin/wine BEL.exe '$i'.BEL_config.'$BEL_CONFIG_FILE_SUFFIX'.txt"'
 echo cd ..
end
