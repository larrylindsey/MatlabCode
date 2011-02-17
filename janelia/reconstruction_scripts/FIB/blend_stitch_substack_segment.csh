#!/bin/csh

echo "Blend stitched segmentation from overlapping sub-stacks into"
echo "one big stack."
if($#argv == 0) then
    echo "usage: $0 output_tif_name"
    exit
endif

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_FILTER_VERSION=$argv[2]

# segmentation parameters
set SUPERPIXEL_SUFFIX=$argv[3]
set SEGMENT_SUFFIX=$SUPERPIXEL_SUFFIX$argv[4]

set BLEND_OUTPUT=$argv[5]

set OUTPUT_SCALE=$argv[6]

# stitching parameters
set STITCH_SUFFIX=".$argv[7]_$argv[$#argv]"

set STACK_SECTION_SEQ="$argv[7-]"

# constants
set EM_ROOT_DIR="/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline"
echo EM_ROOT_DIR=$EM_ROOT_DIR
set CODE_DIR=$EM_ROOT_DIR"/code/segmentation_3D"
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/"$STACK_NAME
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"
set STITCH_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/stitch"

# derived variables
set IMAGE_SUB_STACK_SUFFIX = ".%d.%d."$STACK_FILTER_VERSION
set IMAGE_SUB_STACK = "image_stack"$IMAGE_SUB_STACK_SUFFIX
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"
set STITCH_OUTPUT=$STITCH_DIR/"stitch_seg"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX$STITCH_SUFFIX".raw"


# stuff
echo $CODE_DIR/blend_stitch_substack_segment_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $STITCH_OUTPUT $BLEND_OUTPUT $STACK_SECTION_SEQ
$CODE_DIR/blend_stitch_substack_segment_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $STITCH_OUTPUT $BLEND_OUTPUT $OUTPUT_SCALE $STACK_SECTION_SEQ
