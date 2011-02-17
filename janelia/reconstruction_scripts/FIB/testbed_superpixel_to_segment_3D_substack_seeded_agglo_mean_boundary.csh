#!/bin/csh

#set i=46
#set section_sequence=""
#set substack_depth=30
#while ($i<76) #135 192 #472
#    @ j = $i + $substack_depth - 1
#    set section_sequence="$section_sequence $i $j"
#    @ i = $i + $substack_depth - 1
#end

set section_sequence="1 30 30 59 59 84"

set s = ($section_sequence)
@ n_substack= $#s / 2

set stack_filter="mf5"
set f_threshold_ws_ld=60
set area_threshold_ld=1000
set superpixel_suffix=".ws.T$f_threshold_ws_ld"
set segment_method="ladder"
set segment_suffix_ld=".ld.T"$f_threshold_ws_ld"_L"$area_threshold_ld
set f_threshold_seq_amb="140 160 180 200"
set segment_suffix_amb="$segment_suffix_ld.ambcs.T200.46_104"
set output_suffix_amb="."$stack_filter$superpixel_suffix$segment_suffix_amb

echo ./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "rat_cortex.FIB.013009" $stack_filter $segment_method $superpixel_suffix $segment_suffix_ld $n_substack $section_sequence $f_threshold_seq_amb
./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix $segment_method $segment_suffix_ld $n_substack $section_sequence $f_threshold_seq_amb

echo ./blend_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence
./blend_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence
