CC = gcc
CPP = gcc
EMLIB_DIR = ../.
EMBIN = .
MYERS_DIR = $(EMLIB_DIR)/genelib
INCLUDE_DIR = -I/usr/local/include -I$(MYERS_DIR) -I./ -I$(EMLIB_DIR)/stack_3D_utilities -I$(EMLIB_DIR)/merge_sets -I$(EMLIB_DIR)/hash_functions -I$(EMLIB_DIR)/image_processing_utilities
IS_ENABLED_DEBUG = __NO_DEBUG__
CFLAGS = $(INCLUDE_DIR) -g -O3 -fPIC -Wno-deprecated -D$(IS_ENABLED_DEBUG)
CPPFLAGS = $(INCLUDE_DIR) -g -O3 -fPIC -Wno-deprecated
LIB_DIR = -L$(MYERS_DIR) -L/usr/local/lib -L$(EMLIB_DIR)/merge_sets -L$(EMLIB_DIR)/image_processing_utilities -L./
LIB =  -liput -lmyers -lm -lsegment
LDFLAGS = $(LIB_DIR) $(LIB) -pg

libsegment_modules = segment_agglo_RAG.o agglo_merge_criterion_mean_boundary.o agglo_merge_criterion_mean_boundary_ladder.o 

all: libsegment.a segment_agglo_RAG_mean_boundary_ladder

segment_agglo_RAG_mean_boundary_ladder: segment_agglo_RAG_mean_boundary_ladder.o libsegment.a
	$(CC) $@.o $(LDFLAGS) -lstdc++ -o $(EMBIN)/$@
	rm $@.o

libsegment.a: $(libsegment_modules)
	ar -r libsegment.a $(libsegment_modules)
	rm $(libsegment_modules)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm *.o
