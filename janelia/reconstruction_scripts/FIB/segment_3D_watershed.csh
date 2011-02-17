#!/bin/csh

echo "Compute 3D watershed superpixels"

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_START_SECTION=$argv[2]
set STACK_END_SECTION=$argv[3]
set STACK_FILTER_VERSION=$argv[4]

# segmentation parameters
set F_THRESHOLD=$argv[5]

# constants
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/$STACK_NAME"
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"

# derived variables
set IMAGE_SUB_STACK_SUFFIX = "."$STACK_START_SECTION"."$STACK_END_SECTION"."$STACK_FILTER_VERSION
set IMAGE_SUB_STACK = "image_stack"$IMAGE_SUB_STACK_SUFFIX
set SUPERPIXEL_SUFFIX=".ws.T"$F_THRESHOLD
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"

# stuff
echo $CODE_DIR/watershed_3D_myers $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $F_THRESHOLD $SUPERPIXEL
$CODE_DIR/watershed_3D_myers $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $F_THRESHOLD $SUPERPIXEL
