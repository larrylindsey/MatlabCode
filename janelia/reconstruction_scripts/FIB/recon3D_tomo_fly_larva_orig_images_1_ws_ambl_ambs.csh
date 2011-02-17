#!/bin/csh

setenv EM_ROOT_DIR /groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline
setenv CODE_DIR /groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/code/segmentation_3D

set is_enabled_watershed=1
set is_enabled_agglo_mean_boundary_ladder=1
set is_enabled_agglo_mean_boundary_ladder_blend=0
set is_enabled_combined_agglo_mean_boundary=1
set is_enabled_combined_agglo_mean_boundary_blend=0
set is_enabled_combined_seeded_agglo_mean_boundary=0
set is_enabled_combined_seeded_agglo_mean_boundary_blend=0

set stack_name="tomo_fly_larva_orig_images_1"
#set i=1
#set section_sequence=""
#set substack_depth=10
#while ($i<9)
#    @ j = $i + $substack_depth - 1
#    set section_sequence="$section_sequence $i $j"
#    @ i = $i + $substack_depth - 1
#end

set section_sequence="1 30"

set stack_filter="heq_d4_mf8"
set f_threshold_ws=1
set superpixel_suffix=".ws.T$f_threshold_ws"
set s = ($section_sequence)
set f_threshold_ambl=120
set boundary_length_threshold_ambl=50
set minimum_area_threshold_ambl=1000
set f_threshold_seq_ambcs="130 140 150 160 170 190 200 210 220 230 240 250"
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
        echo ./segment_3D_watershed.csh "$stack_name" $i $j $stack_filter $f_threshold_ws 
        ./segment_3D_watershed.csh "$stack_name" $i $j $stack_filter $f_threshold_ws  
    endif

    if($is_enabled_agglo_mean_boundary_ladder == 1) then
        echo ./superpixel_to_segment_3D_agglo_mean_boundary_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_ambl $f_threshold_ambl $minimum_area_threshold_ambl
        ./superpixel_to_segment_3D_agglo_mean_boundary_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_ambl $f_threshold_ambl $minimum_area_threshold_ambl
    endif

    if($is_enabled_agglo_mean_boundary_ladder_blend == 1) then
        set segment_suffix_ambl=".ambl.T$f_threshold_ambl.L$boundary_length_threshold_ambl.A$minimum_area_threshold_ambl"
        set output_suffix_ambl="."$stack_filter$superpixel_suffix$segment_suffix_ambl".$i.$j"
        echo ./blend_stack_segment.csh "$stack_name" $i $j $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl "~/temp/fib.seg"$output_suffix_ambl".tif"
        ./blend_stack_segment.csh "$stack_name" $i $j $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl "~/temp/fib.seg"$output_suffix_ambl".tif"
    endif
end

if($is_enabled_combined_agglo_mean_boundary == 1) then
    set segment_suffix_ambcs=".ambl.T$f_threshold_ambl.L$boundary_length_threshold_ambl.A$minimum_area_threshold_ambl"
    echo ./superpixel_to_segment_3D_substack_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs $n_substack $section_sequence $f_threshold_seq_ambcs
    ./superpixel_to_segment_3D_substack_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs $n_substack $section_sequence $f_threshold_seq_ambcs
endif

if($is_enabled_combined_seeded_agglo_mean_boundary == 1) then
    set segment_suffix_ambcs=".ambl.T$f_threshold_ambl.L$boundary_length_threshold_ambl.A$minimum_area_threshold_ambl"
    echo ./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs $n_substack $section_sequence $f_threshold_seq_ambcs
    ./superpixel_to_segment_3D_substack_seeded_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs $n_substack $section_sequence $f_threshold_seq_ambcs
endif

if($is_enabled_combined_seeded_agglo_mean_boundary_blend == 1) then
    foreach f ($f_thresholds_ambcs_blend)
        set segment_suffix_ambl=".ambl.T$f_threshold_ambl.L$boundary_length_threshold_ambl.A$minimum_area_threshold_ambl"
        set segment_suffix_ambcs=".ambcs.T$f.$s[1]_$s[$#s]"
        set output_suffix_ambcs=".$stack_filter$superpixel_suffix$segment_suffix_ambl$segment_suffix_ambcs"
        echo ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl "~/temp/fib.seg."$output_suffix_ambcs".tif" 2 $section_sequence
        ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl "~/temp/fib.seg."$output_suffix_ambcs".tif" 1 $section_sequence
    end
endif
