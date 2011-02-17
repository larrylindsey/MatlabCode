#!/bin/csh

echo "Stitch together segmentation from overlapping sub-stacks"
#if($#argv == 0) then
#    echo "usage: $0"
#    exit
#endif

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_FILTER_VERSION=$argv[2]

# segmentation parameters
set SUPERPIXEL_SUFFIX=$argv[3]
set SEGMENT_METHOD=$argv[4]
set SEGMENT_SUFFIX=$SUPERPIXEL_SUFFIX$argv[5]

set STACK_SECTION_SEQ="$argv[6-]"

# stitching parameters
set AREA_OVERLAP_THRESHOLD=3000
set AREA_OVERLAP_NORM1_THRESHOLD=0.15
set AREA_OVERLAP_NORM2_THRESHOLD=0.7
set STITCH_SUFFIX=".$argv[6]_$argv[$#argv]"

# constants
set EM_ROOT_DIR="/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline"
echo EM_ROOT_DIR=$EM_ROOT_DIR
set CODE_DIR=$EM_ROOT_DIR"/code/segmentation_3D"
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/"$STACK_NAME
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"
set SEGMENT_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/$SEGMENT_METHOD"
set STITCH_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/stitch"

# derived variables
set IMAGE_SUB_STACK_SUFFIX = ".%d.%d."$STACK_FILTER_VERSION
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"
set SEGMENT=$SEGMENT_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX".raw"
set STITCH_OUTPUT=$STITCH_DIR/"stitch_seg"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX$STITCH_SUFFIX".raw"

# stuff
echo $CODE_DIR/stitch_substack_segment_overlap_aov0_b $SUPERPIXEL $SEGMENT $AREA_OVERLAP_THRESHOLD $AREA_OVERLAP_NORM1_THRESHOLD $AREA_OVERLAP_NORM2_THRESHOLD $STITCH_OUTPUT $STACK_SECTION_SEQ
$CODE_DIR/stitch_substack_segment_overlap_aov0_b $SUPERPIXEL $SEGMENT $AREA_OVERLAP_THRESHOLD $AREA_OVERLAP_NORM1_THRESHOLD $AREA_OVERLAP_NORM2_THRESHOLD $STITCH_OUTPUT $STACK_SECTION_SEQ


