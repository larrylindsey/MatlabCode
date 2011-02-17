#!/bin/csh

echo "Compute 3D ladder"

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_START_SECTION=$argv[2]
set STACK_END_SECTION=$argv[3]
set STACK_FILTER_VERSION=$argv[4]

# segmentation parameters
set SUPERPIXEL_SUFFIX=$argv[5]
set F_THRESHOLD=$argv[6]
set AREA_THRESHOLD=$argv[7]

# constants
set EM_ROOT_DIR="/groups/chklovskii/home/vitaladevunis/research/em_reconstruction_pipeline"
set CODE_DIR=$EM_ROOT_DIR"/code/segmentation_3D"
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/"$STACK_NAME
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"
set SEGMENT_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/ladder"

# derived variables
set IMAGE_SUB_STACK_SUFFIX = "."$STACK_START_SECTION"."$STACK_END_SECTION"."$STACK_FILTER_VERSION
set IMAGE_SUB_STACK = "image_stack"$IMAGE_SUB_STACK_SUFFIX
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"
set SEGMENT=$SEGMENT_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".ld.T"$F_THRESHOLD"_L"$AREA_THRESHOLD".raw"

# stuff
echo $CODE_DIR/superpixel_2_segment_3D_ladder_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $F_THRESHOLD $AREA_THRESHOLD $SEGMENT
$CODE_DIR/superpixel_2_segment_3D_ladder_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $F_THRESHOLD $AREA_THRESHOLD $SEGMENT
