#!/bin/csh
set BEL_CONFIG_FILE_SUFFIX=$1
if ($?WINE_PATH == 0)  set WINE_PATH=/usr/local/wine/bin

foreach i ( $argv[2-])
 echo cd $i
 echo 'echo "'$WINE_PATH'/wine BEL.exe '$i'.BEL_config.'$BEL_CONFIG_FILE_SUFFIX'.txt" | batch'
 echo cd ..
end
