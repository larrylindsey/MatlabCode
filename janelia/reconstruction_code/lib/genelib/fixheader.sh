#!/bin/bash

#this script is used to fix absolute include path

for file in `ls *.h` 
do 
    sed 's|/Users/myersg/src/||' <$file >$file.new 
    mv $file.new $file
done

for file in `ls *.p` 
do 
    sed 's|/Users/myersg/src/||' <$file >$file.new 
    mv $file.new $file
done