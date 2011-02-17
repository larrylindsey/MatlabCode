#!/bin/csh

if (! -d $argv[1]) then
    echo "creating dir "$argv[1]
    mkdir $argv[1]
endif
