#!/bin/csh

set i=46
set section_sequence=""
set section_start=""
set substack_depth=30
while ($i<76) #135 192 #472
    @ j = $i + $substack_depth - 1
    set section_sequence="$section_sequence $i $j"
    set section_start="$section_start $i"
    @ i = $i + $substack_depth - 1
end

set stack_filter="mf5"
set f_threshold_ws_ld=60
set area_threshold_ld=1000
set superpixel_suffix=".ws.T$f_threshold_ws_ld"
set segment_suffix_ld=".ld.T"$f_threshold_ws_ld"_L"$area_threshold_ld
set output_suffix_ld="."$stack_filter$superpixel_suffix$segment_suffix_ld".46.162"
set f_threshold_seq_amb="70 80 90 100 110 120 130 140 150 160"
set segment_suffix_amb=$segment_suffix_ld".amb.T120"
set output_suffix_amb="."$stack_filter$superpixel_suffix$segment_suffix_amb".46.162"

foreach i ($section_start)
    @ j = $i + $substack_depth - 1

#    echo ./segment_3D_watershed.csh "rat_cortex.FIB.013009" $i $j $stack_filter $f_threshold_ws_ld   
#    ./segment_3D_watershed.csh "rat_cortex.FIB.013009" $i $j $stack_filter $f_threshold_ws_ld   

#    echo ./superpixel_to_segment_3D_ladder.csh "rat_cortex.FIB.013009" $i $j $stack_filter $superpixel_suffix $f_threshold_ws_ld $area_threshold_ld
#    ./superpixel_to_segment_3D_ladder.csh "rat_cortex.FIB.013009" $i $j $stack_filter $superpixel_suffix $f_threshold_ws_ld $area_threshold_ld

#    echo ./blend_stack_segment.csh "rat_cortex.FIB.013009" $i $j $stack_filter $superpixel_suffix "ladder" $segment_suffix_ld "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_ld".tif"
#    ./blend_stack_segment.csh "rat_cortex.FIB.013009" $i $j "$stack_filter" $superpixel_suffix "ladder" $segment_suffix_ld "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_ld".tif"

#    echo ./superpixel_to_segment_3D_agglo_mean_boundary.csh "rat_cortex.FIB.013009" $i $j "$stack_filter" $superpixel_suffix $segment_suffix_ld $f_threshold_seq_amb
#    ./superpixel_to_segment_3D_agglo_mean_boundary.csh "rat_cortex.FIB.013009" $i $j "$stack_filter" $superpixel_suffix $segment_suffix_ld $f_threshold_seq_amb

#    echo ./blend_stack_segment.csh "rat_cortex.FIB.013009" $i $j "$stack_filter" $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_amb".tif"
#    ./blend_stack_segment.csh "rat_cortex.FIB.013009" $i $j "$stack_filter" $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_amb".tif"
end

## For ladder ##
#echo ./stitch_substack_segment_overlap.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "ladder" $segment_suffix_ld $section_sequence
#./stitch_substack_segment_overlap.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "ladder" $segment_suffix_ld $section_sequence

#echo ./blend_stitch_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix $segment_suffix_ld "~/temp/fib.seg."$output_suffix_ld".tif" 2 $section_sequence
#./blend_stitch_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix $segment_suffix_ld "~/temp/fib.seg."$output_suffix_ld".tif" 2 $section_sequence

## For agglomerative mean boundary ##
echo ./stitch_substack_segment_overlap.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb $section_sequence
./stitch_substack_segment_overlap.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb $section_sequence

echo ./blend_stitch_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence
./blend_stitch_substack_segment.csh "rat_cortex.FIB.013009" $stack_filter $superpixel_suffix $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence

