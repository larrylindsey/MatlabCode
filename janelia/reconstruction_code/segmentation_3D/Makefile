CC = gcc
CPP = gcc
EMLIB_DIR = ../lib
EMBIN = .
MYERS_DIR = $(EMLIB_DIR)/genelib
INCLUDE_DIR = -I/usr/local/include -I$(MYERS_DIR) -I./ -I$(EMLIB_DIR)/stack_3D_utilities -I$(EMLIB_DIR)/merge_sets -I$(EMLIB_DIR)/image_processing_utilities
CFLAGS = $(INCLUDE_DIR) -g -fPIC -O3 -Wno-deprecated
CPPFLAGS = $(INCLUDE_DIR) -g -fPIC -O3 -Wno-deprecated
LIB_DIR = -L$(MYERS_DIR) -L/usr/local/lib -L$(EMLIB_DIR)/merge_sets -L$(EMLIB_DIR)/image_processing_utilities
LIB =  -liput -lmyers -lm
LDFLAGS = $(LIB_DIR) $(LIB) -pg

all: blend_stitch_substack_segment blend_stack_segment watershed_3D_myers superpixel_2_segment_3D_ladder superpixel_2_segment_3D_agglo_mean_boundary superpixel_2_segment_3D_substack_agglo_mean_boundary superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary stitch_substack_segment_overlap_aov0 libsegment_3D_utilities.a

watershed_3D_myers: watershed_3D.o $(EMLIB_DIR)/image_processing_utilities/libiput.a $(MYERS_DIR)/libmyers.a
	$(CC) watershed_3D.o $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@
	@rm watershed_3D.o

superpixel_2_segment_3D_ladder: superpixel_2_segment_3D_ladder.o $(EMLIB_DIR)/merge_sets/libmerge_sets.a $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lmerge_sets -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

superpixel_2_segment_3D_agglo_mean_boundary: superpixel_2_segment_3D_agglo_mean_boundary.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

superpixel_2_segment_3D_agglo_mean_boundary_ladder: superpixel_2_segment_3D_agglo_mean_boundary_ladder.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

superpixel_2_segment_3D_substack_agglo_mean_boundary: superpixel_2_segment_3D_substack_agglo_mean_boundary.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary: superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

blend_stack_segment: blend_stack_segment.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

blend_stitch_substack_segment: blend_stitch_substack_segment.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@_b
	@rm $@.o

stitch_substack_segment_overlap_aov0: stitch_substack_segment_overlap_aov0.o $(MYERS_DIR)/libmyers.a $(EMLIB_DIR)/image_processing_utilities/libiput.a
	$(CC) $@.o  $(LDFLAGS) -lstdc++ -lmerge_sets -o $(EMBIN)/$@_b
	@rm $@.o

libsegment_3D_utilities_modules = remove_merged_boundaries.o draw_segment_boundaries.o

libsegment_3D_utilities.a: $(libsegment_3D_utilities_modules)
	ar -r libsegment_3D_utilities.a $(libsegment_3D_utilities_modules)
	rm $(libsegment_3D_utilities_modules)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
