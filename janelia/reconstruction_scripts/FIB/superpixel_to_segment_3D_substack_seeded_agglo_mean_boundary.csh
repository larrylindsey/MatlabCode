#!/bin/csh

echo "Compute 3D agglo. mean boundary. Combine multiple substacks."

# stack specific parameters
set STACK_NAME=$argv[1]
set STACK_FILTER_VERSION=$argv[2]

# segmentation parameters
set SUPERPIXEL_SUFFIX=$argv[3]
set SEGMENT_METHOD=$argv[4]
set SEGMENT_SUFFIX=$SUPERPIXEL_SUFFIX$argv[5]
set ARGC_N_SUBSTACK=6
set N_SUBSTACK=$argv[$ARGC_N_SUBSTACK]
echo "N_SUBSTACK: $N_SUBSTACK"
@ ARGC_SUBSTACK_SEQ_START = $ARGC_N_SUBSTACK + 1
@ ARGC_SUBSTACK_SEQ_END = $ARGC_N_SUBSTACK + $N_SUBSTACK + $N_SUBSTACK
set SUBSTACK_SEQ="$argv[$ARGC_SUBSTACK_SEQ_START-$ARGC_SUBSTACK_SEQ_END]"
echo "SUBSTACK_SEQ: $SUBSTACK_SEQ"
set STACK_START_SECTION=$argv[$ARGC_SUBSTACK_SEQ_START]
set STACK_END_SECTION=$argv[$ARGC_SUBSTACK_SEQ_END]
echo "STACK_START_SECTION: $STACK_START_SECTION, STACK_END_SECTION: $STACK_END_SECTION"
@ ARGC_F_THRESHOLD_SEQ = 1 + $ARGC_SUBSTACK_SEQ_END
set F_THRESHOLD_SEQ="$argv[$ARGC_F_THRESHOLD_SEQ-]"

# constants
set RECONSTRUCTION_DIR=$EM_ROOT_DIR"/reconstructions/"$STACK_NAME
set SUB_STACK_DIR=$RECONSTRUCTION_DIR"/image_stacks"
set SUPERPIXEL_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/watershed"
set SEGMENT_INPUT_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/$SEGMENT_METHOD"
set SEGMENT_OUTPUT_DIR=$RECONSTRUCTION_DIR"/3D_segmentation_results/agglo_mean_boundary"

# derived variables
set IMAGE_SUB_STACK_SUFFIX = ".%u.%u."$STACK_FILTER_VERSION
set IMAGE_SUB_STACK = "image_stack"$IMAGE_SUB_STACK_SUFFIX
set SUPERPIXEL=$SUPERPIXEL_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SUPERPIXEL_SUFFIX".raw"
set SEGMENT_INPUT=$SEGMENT_INPUT_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX".raw"
set SEGMENT_OUTPUT=$SEGMENT_OUTPUT_DIR/"seg_stack"$IMAGE_SUB_STACK_SUFFIX$SEGMENT_SUFFIX".ambcs.T%u.$argv[$ARGC_SUBSTACK_SEQ_START]_$argv[$ARGC_SUBSTACK_SEQ_END].raw"

# stuff
echo $CODE_DIR/superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $SEGMENT_INPUT $SEGMENT_OUTPUT $N_SUBSTACK $SUBSTACK_SEQ $F_THRESHOLD_SEQ
$CODE_DIR/superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary_b $SUB_STACK_DIR/$IMAGE_SUB_STACK".tif" $SUPERPIXEL $SEGMENT_INPUT $SEGMENT_OUTPUT $N_SUBSTACK $SUBSTACK_SEQ $F_THRESHOLD_SEQ

