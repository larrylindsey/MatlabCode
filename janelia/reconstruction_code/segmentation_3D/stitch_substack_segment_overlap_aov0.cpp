// Combine two overlapping segment maps by remapping their labels to a common
// range.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	03202009	init. code
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <utility>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>

#include <image_lib.h>

#include <file_input_output.h>
#include <merge_sets.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

struct ltstru{
  bool operator()(const unsigned int l1, unsigned int l2) const
    {
      return l1 < l2;
    }
};
struct ltstrlu{
  bool operator()(const unsigned long int l1, unsigned long int l2) const
    {
      return l1 < l2;
    }
};

#define ARGC_SUPERPIXEL_NAME_FORMAT 1
#define ARGC_SEGMENT_NAME_FORMAT 2
#define ARGC_AREA_OVERLAP_THRESHOLD 3
#define ARGC_AREA_OVERLAP_NORM1_THRESHOLD 4
#define ARGC_AREA_OVERLAP_NORM2_THRESHOLD 5
#define ARGC_OUTPUT_NAME_FORMAT 6
#define ARGC_FIRST_SECTION 7
int main(int argc, char * argv[]){
	if(argc==0){
		printf("Usage stitch_substack_segment_overlap_aov0_b superpixel_name_format segment_name_format area_overlap_threshold area_overlap_norm1_threshold area_overlap_norm2_threshold output_name_format section_sequence\n");
		printf("input params:\n");
		printf("\t1. Format of superpixel file names.\n");
    printf("\t2. Format of segment file names.\n");
		printf("\t3. Rule1 area of overlap threshold (pixels)\n");
		printf("\t4. Rule1 normalized area of overlap threshold\n");
		printf("\t5. Rule2 normalized area of overlap threshold\n");
		printf("\t6. output label file name format\n");
    printf("\t7. sequence of sub-section start and end.\n");
		return 1;
	}

  printf("\nSTART: combine_segment_maps_3D_aov0_m2m_b\n");
  printf("superpixel file name format: %s\n",
         argv[ARGC_SUPERPIXEL_NAME_FORMAT]);
  printf("segment file name format: %s\n",
         argv[ARGC_SEGMENT_NAME_FORMAT]);
  int n_substack = 0;
  {
    printf("Sub-stack begin and end:\n");
    int i;
    for(i=ARGC_FIRST_SECTION; i<argc; i+=2){
      printf("%d %d\n", atoi(argv[i]), atoi(argv[i+1]));
      n_substack++;
    }
  }

  int pair_id;
  char * file_name = (char *) new char
    [MAX(MAX(strlen(argv[ARGC_SUPERPIXEL_NAME_FORMAT]),
             strlen(argv[ARGC_SEGMENT_NAME_FORMAT])),
         strlen(argv[ARGC_OUTPUT_NAME_FORMAT]))];
                     
  unsigned int * superpixel_map_1 = NULL;
  int stack_width_1;
  int stack_height_1;
  int stack_depth_1;
  int stack_n_pixel_1;
  int stack_n_pixel_plane_1;
  unsigned int max_superpixel_label_1 = 0;
  unsigned int * label_mapping_1;
  unsigned int max_label_mapping_1=0;
  Merge_Sets label_sets(0, NULL);
  unsigned int * label_sets_offsets = (unsigned int *) new
    unsigned int[n_substack+1];
  label_sets_offsets[0] = 0;
  for(pair_id=0; pair_id<n_substack-1; pair_id++){
    printf("pair: %d %d\n", pair_id, pair_id+1);
    int argc_pair_1 = ARGC_FIRST_SECTION + pair_id*2;
    int argc_pair_2 = ARGC_FIRST_SECTION + pair_id*2 + 2;
    int start_section_1 = atoi(argv[argc_pair_1]);
    int end_section_1 = atoi(argv[argc_pair_1+1]);    
    int start_section_2 = atoi(argv[argc_pair_2]);
    int end_section_2 = atoi(argv[argc_pair_2+1]);
    
    if(superpixel_map_1==NULL){
      printf("Reading in the superpixel stack %d\n", pair_id);
      {
        sprintf(file_name, argv[ARGC_SUPERPIXEL_NAME_FORMAT],
                start_section_1, end_section_1);
        iput::io::Type_Element type_element;
        int n_dimension, * dimensions;
        void * read_data = iput::io::fread_raw_array(file_name,
                                                     type_element, n_dimension,
                                                     dimensions);
        if(iput::io::is_error_type_element(type_element) || read_data==NULL){
          printf("Could not open superpixel map\n");
          return -1;
        }
        if(n_dimension!=3){
          printf("Not a 3D segment stack\n");
          return -1;
        }
        stack_width_1 = dimensions[0];
        stack_height_1 = dimensions[1];
        stack_depth_1 = dimensions[2];
        stack_n_pixel_1 = stack_width_1*stack_height_1*stack_depth_1;
        stack_n_pixel_plane_1 = stack_width_1*stack_height_1;
        superpixel_map_1 = (unsigned int *) read_data;
        delete [] dimensions;
      }
      printf("stack_width_1=%d, stack_height_1=%d, stack_depth_1=%d, stack_n_pixel_1=%d, stack_n_pixel_plane_1=%d\n",
             stack_width_1, stack_height_1, stack_depth_1, stack_n_pixel_1,
             stack_n_pixel_plane_1);

      {
        int i;
        unsigned int * l = superpixel_map_1;
        for(i=0; i<stack_n_pixel_1; i++, l++)
          max_superpixel_label_1 = MAX(max_superpixel_label_1, *l);
      }
  
      printf("Reading in the segment mapping %d\n", pair_id);
      {
        sprintf(file_name, argv[ARGC_SEGMENT_NAME_FORMAT],
                start_section_1, end_section_1);
        iput::io::Type_Element type_element;
        int n_dimension, * dimensions;
        void * read_data = iput::io::fread_raw_array(file_name,
                                                     type_element, n_dimension,
                                                     dimensions);
        if(iput::io::is_error_type_element(type_element) || read_data==NULL){
          printf("\nERROR: Could not open label mapping\n");
          return -1;
        }
        if(n_dimension!=1){
          printf("\nERROR: Not a 1D label mapping\n");
          return -1;
        }
        if(dimensions[0]!=max_superpixel_label_1+1){
          printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
          return -1;
        }
        label_mapping_1 = (unsigned int *) read_data;
        for(int i=0; i<max_superpixel_label_1; i++)
          max_label_mapping_1 = MAX(max_label_mapping_1, label_mapping_1[i]);
        delete [] dimensions;
      }
      
      label_sets.add_new_sets(max_label_mapping_1);
    }
    printf("max_superpixel_label_1: %d\n", max_superpixel_label_1);
    printf("max_label_mapping_1: %d\n", max_label_mapping_1);
    printf("no. of sets: %d\n", label_sets.adam.size());
    
    printf("Reading in the superpixel stack %d\n", pair_id);
    int stack_width_2;
    int stack_height_2;
    int stack_depth_2;
    int stack_n_pixel_2;
    int stack_n_pixel_plane_2;
    unsigned int * superpixel_map_2;
    {
      sprintf(file_name, argv[ARGC_SUPERPIXEL_NAME_FORMAT],
              start_section_2, end_section_2);
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      void * read_data = iput::io::fread_raw_array(file_name,
                                                   type_element, n_dimension,
                                                   dimensions);
      if(iput::io::is_error_type_element(type_element) || read_data==NULL){
        printf("Could not open superpixel map\n");
        return -1;
      }
      if(n_dimension!=3){
        printf("Not a 3D segment stack\n");
        return -1;
      }
      stack_width_2 = dimensions[0];
      stack_height_2 = dimensions[1];
      stack_depth_2 = dimensions[2];
      stack_n_pixel_2 = stack_width_2*stack_height_2*stack_depth_2;
      stack_n_pixel_plane_2 = stack_width_2*stack_height_2;
      superpixel_map_2 = (unsigned int *) read_data;
      delete [] dimensions;
    }
    printf("stack_width_2=%d, stack_height_2=%d, stack_depth_2=%d, stack_n_pixel_2=%d, stack_n_pixel_plane_2=%d\n",
           stack_width_2, stack_height_2, stack_depth_2, stack_n_pixel_2, stack_n_pixel_plane_2);
    if(stack_width_1!=stack_width_2 || stack_height_1!=stack_height_2){
      printf("\nERROR: Stack planes do not match.\n");
      return -1;
    }

    unsigned int max_superpixel_label_2 = 0;
    {
      int i;
      unsigned int * l = superpixel_map_2;
      for(i=0; i<stack_n_pixel_2; i++, l++)
        max_superpixel_label_2 = MAX(max_superpixel_label_2, *l);
    }
  
    printf("Reading in the segment mapping %d\n", pair_id);
    unsigned int * label_mapping_2;
    unsigned int max_label_mapping_2=0;
    {
      sprintf(file_name, argv[ARGC_SEGMENT_NAME_FORMAT],
              start_section_2, end_section_2);
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      void * read_data = iput::io::fread_raw_array(file_name,
                                                   type_element, n_dimension,
                                                   dimensions);
      if(iput::io::is_error_type_element(type_element) || read_data==NULL){
        printf("\nERROR: Could not open label mapping\n");
        return -1;
      }
      if(n_dimension!=1){
        printf("\nERROR: Not a 1D label mapping\n");
        return -1;
      }
      if(dimensions[0]!=max_superpixel_label_2+1){
        printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
        return -1;
      }
      label_mapping_2 = (unsigned int *) read_data;
      {
        unsigned int m_l = 0;
        for(int i=0; i<max_superpixel_label_2; i++)
          m_l = MAX(m_l, label_mapping_2[i]);
        printf("initial max_label_mapping_2: %d\n", m_l);
      }
      for(int i=0; i<max_superpixel_label_2; i++){
        label_mapping_2[i] += max_label_mapping_1 + 1;
        max_label_mapping_2 = MAX(max_label_mapping_2, label_mapping_2[i]);
      }
      delete [] dimensions;
    }
    
    label_sets_offsets[pair_id+1] = max_label_mapping_1+1;
    // create new sets for each new segment
    label_sets.add_new_sets(max_label_mapping_2 - max_label_mapping_1);
    // merge the 0 label segment with the 1st sub-stack's 0 label
    label_sets.merge(0, max_label_mapping_1 + 1);
    
    printf("max_superpixel_label_2: %d\n", max_superpixel_label_2);
    printf("max_label_mapping_2: %d\n", max_label_mapping_2);
    printf("no. of sets: %d\n", label_sets.adam.size());

    
    //////////////////////////////////
    // Get the area of overlap
    //////////////////////////////////
    printf("Computing area of labels and overlap\n");
    int n_section_overlap = end_section_1 - start_section_2 + 1;
    printf("n_section_overlap: %d\n", n_section_overlap);
    if(n_section_overlap<=0)
      return -1;
    
    std::map<const unsigned int, int, ltstru> label_area_1, label_area_2;
    std::map<const long unsigned int, int, ltstrlu> label_pair_area_overlap;
    unsigned int max_label_overlap_1=0, max_label_overlap_2=0;
    {
      unsigned int *l1 = superpixel_map_1 + stack_n_pixel_plane_1*(stack_depth_1-n_section_overlap);
      unsigned int *l2 = superpixel_map_2;
      unsigned int s1, s2;
      int i;
      long unsigned int l_p, l;
      for(i=0; i<n_section_overlap*stack_n_pixel_plane_1; i++, l1++, l2++){
        s1 = label_mapping_1[*l1];
        s2 = label_mapping_2[*l2];
        l_p = s1;
        l = s2;
        l_p |= l << 32;
        label_pair_area_overlap[l_p] = label_pair_area_overlap[l_p] + 1;
        label_area_1[s1]++;
        label_area_2[s2]++;
      }

#ifdef DEBUG_MODE
      // Debug
      i=0;
      std::map<const long unsigned int, int, ltstrlu>::iterator pair_it;
      for(pair_it=label_pair_area_overlap.begin(); pair_it!=label_pair_area_overlap.end(); pair_it++, i++){
        if(i>1000)
          break;
        long unsigned int l_p = (*pair_it).first;
        int n_pair = (*pair_it).second;
        unsigned int s1, s2;
        s1 = l_p & 0x00000000ffffffff;
        s2 = l_p >> 32;
        printf("l1:%u l2:%u n_pair:%d\n", s1, s2, n_pair);
      }
#endif
    }

    printf("Get sets of overlapping labels that should be linked together.\n");
    printf("Method: Merge sets library\n");
    {
      int area_overlap_threshold = atoi(argv[ARGC_AREA_OVERLAP_THRESHOLD]);
      double area_overlap_norm1_threshold =
        atof(argv[ARGC_AREA_OVERLAP_NORM1_THRESHOLD]);
      double area_overlap_norm2_threshold =
        atof(argv[ARGC_AREA_OVERLAP_NORM2_THRESHOLD]);
      printf("area_overlap_threshold: %d\n", area_overlap_threshold);
      printf("area_overlap_norm1_threshold: %f\n",
             area_overlap_norm1_threshold);
      printf("area_overlap_norm2_threshold: %f\n",
             area_overlap_norm2_threshold);
      std::map<const long unsigned int, int, ltstrlu>::iterator pair_it;
      unsigned long int l_p, l;
      for(pair_it=label_pair_area_overlap.begin();
          pair_it!=label_pair_area_overlap.end(); pair_it++){
        long unsigned int l_p = (*pair_it).first;
        double area_overlap = (double) (*pair_it).second;
        unsigned int s1, s2;
        s1 = l_p & 0x00000000ffffffff;
        s2 = l_p >> 32;
        if(s1==0 || s2==0)
          continue;
        double area_1 = (double) label_area_1[s1];
        double area_2 = (double) label_area_2[s2];
#ifdef DEBUG_MODE
        printf("area_1:%g area_2:%g, overlap:%g\n", area_1, area_2,
               area_overlap);
#endif
        bool rule_1 = area_overlap>area_overlap_threshold
          && area_overlap/MAX(area_1, area_2)>area_overlap_norm1_threshold;
        bool rule_2 = area_overlap/MIN(area_1, area_2)>
          area_overlap_norm2_threshold;
#ifdef DEBUG_MODE
        printf("rule_1:%d rule_2:%d\n", rule_1, rule_2);
#endif
        bool rule = rule_1 || rule_2;
        if(rule){
#ifdef DEBUG_MODE
          printf("s1:%d, s2:%d\n", s1, s2);
          printf("a1:%d, a2:%d\n", label_sets.adam[s1], label_sets.adam[s2]);
#endif
          if(s1<s2)
            label_sets.merge(s1, s2);
          else
            label_sets.merge(s2, s1);
#ifdef DEBUG_MODE
          printf("a1:%d, a2:%d\n", label_sets.adam[s1], label_sets.adam[s2]);
#endif
        }
      }
    }

    delete [] superpixel_map_1;
    delete [] label_mapping_1;
    superpixel_map_1 = superpixel_map_2;
    label_mapping_1 = label_mapping_2;
    stack_depth_1 = stack_depth_2;
    stack_n_pixel_1 = stack_n_pixel_2;
    max_superpixel_label_1 = max_superpixel_label_2;
    max_label_mapping_1 = max_label_mapping_2;
    printf("--------------------------------------------------------\n");
  }
  label_sets_offsets[n_substack] = max_label_mapping_1+1;
  delete [] superpixel_map_1;
  delete [] label_mapping_1;

  label_sets.update_adams();

  printf("Saving merged segment labels ...\n");
  int substack_id;
  for(substack_id=0; substack_id<n_substack; substack_id++){
    printf("substack_id:%d\n", substack_id);
    int argc_substack = ARGC_FIRST_SECTION + substack_id*2;
    unsigned int * label_mapping;
    unsigned int max_superpixel_id;
    {
      sprintf(file_name, argv[ARGC_SEGMENT_NAME_FORMAT],
              atoi(argv[argc_substack]), atoi(argv[argc_substack+1]));
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      void * read_data = iput::io::fread_raw_array(file_name,
                                                   type_element, n_dimension,
                                                   dimensions);
      if(iput::io::is_error_type_element(type_element) || read_data==NULL){
        printf("\nERROR: Could not open label mapping\n");
        return -1;
      }
      if(n_dimension!=1){
        printf("\nERROR: Not a 1D label mapping\n");
        return -1;
      }
      label_mapping = (unsigned int *) read_data;
      max_superpixel_id = dimensions[0]-1;
      delete [] dimensions;
    }
    {
      int size_element = sizeof(unsigned int);
      iput::io::Type_Element type_element =
        iput::io::get_type_element(iput::io::type_id_unsigned_int(),
                                   size_element);
      int n_dimension = 1;
      printf("max_superpixel_id: %d\n", max_superpixel_id);
      int dimensions[1] = {max_superpixel_id+1};
      for(int i=0; i<=max_superpixel_id; i++)
        label_mapping[i] = label_sets.adam[label_mapping[i]+
                                           label_sets_offsets[substack_id]];
      int err;
      {
        sprintf(file_name, argv[ARGC_OUTPUT_NAME_FORMAT],
                atoi(argv[argc_substack]), atoi(argv[argc_substack+1]));
        err = iput::io::fwrite_raw_array(file_name, type_element,
                                         n_dimension, dimensions,
                                         (void *) label_mapping);
      }
      if(err!=0){
        printf("ERROR superpixel_2_segment_3D_ladder_b: could not write label mapping [%d]\n", err);
        delete [] label_mapping;
        return -2;
      }
    }
    delete [] label_mapping;
  }
  printf("done.\n");
  printf("\nSTOP: combine_segment_maps_3D_aov0_m2m_b\n");
  return 0;
}
