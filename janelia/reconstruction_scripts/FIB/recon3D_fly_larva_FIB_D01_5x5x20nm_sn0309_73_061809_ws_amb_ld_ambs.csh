#!/bin/csh

set is_enabled_watershed=0
set is_enabled_agglo_mean_boundary=0
set is_enabled_agglo_mean_boundary_blend=1
set is_enabled_ladder=0
set is_enabled_ladder_blend=0
set is_enabled_combined_agglo_mean_boundary=0
set is_enabled_combined_agglo_mean_boundary_blend=0
set is_enabled_combined_seeded_agglo_mean_boundary=0
set is_enabled_combined_seeded_agglo_mean_boundary_blend=0

set stack_name="fly_larva.FIB.D01_5x5x20nm_sn0309_73.061809"
#set i=1
#set section_sequence=""
#set substack_depth=10
#while ($i<9)
#    @ j = $i + $substack_depth - 1
#    set section_sequence="$section_sequence $i $j"
#    @ i = $i + $substack_depth - 1
#end
#set section_sequence="1 30 30 59 59 84"
set section_sequence="1 30"

set stack_filter="mf5_heq"
set f_threshold_ws_ld=1
set area_threshold_ld=1000
set superpixel_suffix=".ws.T$f_threshold_ws_ld"
set s = ($section_sequence)
set f_threshold_seq_amb="65 70 75 80 85 90 95 100 110 120 130 140"
set boundary_length_threshold_amb=50
set f_thresholds_amb_blend="130 140"
set f_threshold_amb=120
set f_threshold_seq_ambcs="100 110 120 130 140 150 160 170 180 190 200"
set f_thresholds_ambcs_blend="140 160 180 200"

@ n_substack= $#s / 2
set substack_id=0
while($substack_id<$n_substack)
    @ start_id = $substack_id * 2 + 1
    @ end_id = $substack_id * 2 + 2
    @ substack_id = $substack_id + 1
    set i=$s[$start_id]
    set j=$s[$end_id]

    if($is_enabled_watershed == 1) then
        echo ./segment_3D_watershed.csh "$stack_name" $i $j $stack_filter $f_threshold_ws_ld   
        ./segment_3D_watershed.csh "$stack_name" $i $j $stack_filter $f_threshold_ws_ld   
    endif

    if($is_enabled_agglo_mean_boundary == 1) then
        echo ./superpixel_to_segment_3D_agglo_mean_boundary.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_amb $f_threshold_seq_amb
        ./superpixel_to_segment_3D_agglo_mean_boundary.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_amb $f_threshold_seq_amb
    endif

    if($is_enabled_agglo_mean_boundary_blend == 1) then
        set segment_suffix_amb=".amb.T$f_threshold_amb.L$boundary_length_threshold_amb"
        set output_suffix_amb="."$stack_filter$superpixel_suffix$segment_suffix_amb".$i.$j"
        echo ./blend_stack_segment.csh "$stack_name" $i $j $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg"$output_suffix_amb".tif"
        ./blend_stack_segment.csh "$stack_name" $i $j $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg"$output_suffix_amb".tif"
    endif

    if($is_enabled_ladder == 1) then
        echo ./superpixel_to_segment_3D_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $f_threshold_ws_ld $area_threshold_ld
        ./superpixel_to_segment_3D_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $f_threshold_ws_ld $area_threshold_ld
    endif

    if($is_enabled_ladder_blend == 1) then
        echo ./blend_stack_segment.csh "$stack_name" $i $j $stack_filter $superpixel_suffix "ladder" $segment_suffix_ld "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_ld".tif"
        ./blend_stack_segment.csh "$stack_name" $i $j "$stack_filter" $superpixel_suffix "ladder" $segment_suffix_ld "~/temp/test."$i"."$j$superpixel_suffix$segment_suffix_ld".tif"
    endif
end

if($is_enabled_combined_agglo_mean_boundary == 1) then
    echo ./superpixel_to_segment_3D_substack_agglo_mean_boundary.csh "$stack_name" $stack_filter "ladder" $superpixel_suffix $segment_suffix_ld $n_substack $section_sequence $f_threshold_seq_amb
    ./superpixel_to_segment_3D_substack_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "ladder" $segment_suffix_ld $n_substack $section_sequence $f_threshold_seq_amb
endif

if($is_enabled_combined_agglo_mean_boundary_blend == 1) then
    foreach f ($f_thresholds_amb_blend)
        set segment_suffix_amb="$segment_suffix_ld.ambc.T$f.$s[1]_$s[$#s]"
        set output_suffix_amb="."$stack_filter$superpixel_suffix$segment_suffix_amb
        echo ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence
        ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb "~/temp/fib.seg."$output_suffix_amb".tif" 2 $section_sequence
    end
endif

if($is_enabled_combined_seeded_agglo_mean_boundary == 1) then
    set segment_suffix_amb="$segment_suffix_ld.amb.T$f_threshold_amb"
    echo ./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb $n_substack $section_sequence $f_threshold_seq_ambcs
    ./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_amb $n_substack $section_sequence $f_threshold_seq_ambcs
endif

if($is_enabled_combined_seeded_agglo_mean_boundary_blend == 1) then
    foreach f ($f_thresholds_ambcs_blend)
        set segment_suffix_amb="$segment_suffix_ld.amb.T$f_threshold_amb"
        set segment_suffix_ambcs="$segment_suffix_amb.ambcs.T$f.$s[1]_$s[$#s]"
        set output_suffix_ambcs="."$stack_filter$superpixel_suffix$segment_suffix_ambcs
        echo ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs "~/temp/fib.seg."$output_suffix_ambcs".tif" 2 $section_sequence
        ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs "~/temp/fib.seg."$output_suffix_ambcs".tif" 1 $section_sequence
    end
endif
