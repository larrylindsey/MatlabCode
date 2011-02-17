#!/bin/csh

echo "Blend segmentation labels with grayscale stack"
if($#argv == 0) then
    echo "usage: $0 blend_output_name"
    exit
endif

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_START_SECTION=$argv[2]
set STACK_END_SECTION=$argv[3]
set STACK_FILTER_VERSION=$argv[4]

# segmentation parameters
set SUPERPIXEL_SUFFIX=$argv[5]
set SEGMENT_METHOD=$argv[6]
set SEGMENT_SUFFIX=$SUPERPIXEL_SUFFIX$argv[7]

set BLEND_OUTPUT=$argv[8]

# constants
echo EM_ROOT_DIR=$EM_ROOT_DIR
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/"$STACK_NAME
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"
set SEGMENT_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/"$SEGMENT_METHOD

# derived variables
set IMAGE_SUB_STACK_SUFFIX = "."$STACK_START_SECTION"."$STACK_END_SECTION"."$STACK_FILTER_VERSION
set IMAGE_SUB_STACK = "image_stack"$IMAGE_SUB_STACK_SUFFIX
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"
set SEGMENT=$SEGMENT_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX".raw"

# stuff
echo $CODE_DIR/blend_stack_segment_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $SEGMENT $BLEND_OUTPUT
$CODE_DIR/blend_stack_segment_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $SEGMENT $BLEND_OUTPUT
