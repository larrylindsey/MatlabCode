#!/bin/csh

#setenv EM_ROOT_DIR /groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline
#setenv CODE_DIR /groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline/code/segmentation_3D
setenv EM_ROOT_DIR /groups/chklovskii/home/nuneziglesiasj/em_reconstructions
setenv CODE_DIR /groups/chklovskii/home/nuneziglesiasj/em_reconstructions/code/segmentation_3D

set is_enabled_watershed=1
set is_enabled_agglo_mean_boundary_ladder=1
set is_enabled_agglo_mean_boundary_ladder_blend=0
set is_enabled_combined_agglo_mean_boundary=1
set is_enabled_combined_agglo_mean_boundary_blend=0
set is_enabled_combined_seeded_agglo_mean_boundary=0
set is_enabled_combined_seeded_agglo_mean_boundary_blend=0

set stack_name="spams_denoising"
#set i=1
#set section_sequence=""
#set substack_depth=10
#while ($i<9)
#    @ j = $i + $substack_depth - 1
#    set section_sequence="$section_sequence $i $j"
#    @ i = $i + $substack_depth - 1
#end
set section_sequence="0 29"
#set section_sequence="0 29 30 59 60 89 90 119 120 149 150 179 180 209 210 239 240 269 270 299 300 329"
foreach stack_filter ( dict4x4x4.iter112589 dict12x12x12.iter500 dict10x10x10.iter200 dict5x5x5.iter30892 dict8x8x1.iter198371 dict6x6x6.iter10243 dict11x11x1.iter62001 dict7x7x7.iter4476 dict15x15x1.iter19107 dict8x8x8.iter1595 dict18x18x1.iter9348 dict9x9x9.iter615 dict19x19x1.iter7400 dict10x10x10.iter242 dict23x23x1.iter2877 dict11x11x11.iter124 dict27x27x1.iter1042 dict32x32x1.iter397 dict12x12x12.iter70 dict36x36x1.iter221 dict4x4x4.iter226302 dict42x42x1.iter80 dict5x5x5.iter61883 dict8x8x1.iter380598 dict6x6x6.iter21068 dict11x11x1.iter127406 dict7x7x7.iter8927 dict15x15x1.iter38415 dict8x8x8.iter3170 dict18x18x1.iter18438 dict9x9x9.iter1177 dict19x19x1.iter14735 dict23x23x1.iter5746 dict10x10x10.iter457 dict27x27x1.iter2112 dict11x11x11.iter234 dict32x32x1.iter781 dict12x12x12.iter133 dict36x36x1.iter436 dict42x42x1.iter159 dict4x4x4.iter455670 dict8x8x1.iter744407 dict5x5x5.iter124319 dict11x11x1.iter252563 dict6x6x6.iter41993 dict15x15x1.iter77290 dict7x7x7.iter17386 dict18x18x1.iter37525 dict8x8x8.iter6540 dict19x19x1.iter29560 dict9x9x9.iter2326 dict23x23x1.iter11651 dict10x10x10.iter875 )
#set stack_filter="dict12x12x12.iter500"
    echo "now segmenting stack filter: " $stack_filter
    set f_threshold_ws=1
    set superpixel_suffix=".ws.T$f_threshold_ws"
    set s = ($section_sequence)
    set f_threshold_ambl=90
    set boundary_length_threshold_ambl=50
    set minimum_area_threshold_ambl=1000
    set f_threshold_seq_ambcs=`seq 70 3 230`
    set f_thresholds_ambcs_blend="150 160 170 180 190 200"

    if (-d $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results") then
        echo ""
    else
        echo $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results"
        mkdir $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results"
    endif

    if (-d $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/watershed") then
        echo ""
    else
        echo $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/watershed"
        mkdir $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/watershed"
    endif

    if (-d $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/agglo_mean_boundary") then
        echo ""
    else
        echo $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/agglo_mean_boundary"
        mkdir $EM_ROOT_DIR/reconstructions/$stack_name"/3D_segmentation_results/agglo_mean_boundary"
    endif


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
            ./segment_3D_watershed.csh "$stack_name" $i $j $stack_filter $f_threshold_ws  > /dev/null
        endif

        if($is_enabled_agglo_mean_boundary_ladder == 1) then
            echo ./superpixel_to_segment_3D_agglo_mean_boundary_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_ambl $f_threshold_ambl $minimum_area_threshold_ambl
            ./superpixel_to_segment_3D_agglo_mean_boundary_ladder.csh "$stack_name" $i $j $stack_filter $superpixel_suffix $boundary_length_threshold_ambl $f_threshold_ambl $minimum_area_threshold_ambl > /dev/null
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
        ./superpixel_to_segment_3D_substack_agglo_mean_boundary.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambcs $n_substack $section_sequence $f_threshold_seq_ambcs > /dev/null
    endif

    if($is_enabled_combined_agglo_mean_boundary_blend == 1) then
        foreach f ($f_thresholds_ambcs_blend)
            set segment_suffix_ambl=".ambl.T$f_threshold_ambl.L$boundary_length_threshold_ambl.A$minimum_area_threshold_ambl"
            set segment_suffix_ambc=".ambc.T$f.$s[1]_$s[$#s]"
            set output_suffix_ambc=".$stack_filter$superpixel_suffix$segment_suffix_ambl$segment_suffix_ambc"
            echo ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl$segment_suffix_ambc "~/temp/fib.seg."$output_suffix_ambc".tif" 2 $section_sequence
            ./blend_substack_segment.csh "$stack_name" $stack_filter $superpixel_suffix "agglo_mean_boundary" $segment_suffix_ambl$segment_suffix_ambc "~/temp/fib.seg."$output_suffix_ambc".tif" 1 $section_sequence
        end
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

end
