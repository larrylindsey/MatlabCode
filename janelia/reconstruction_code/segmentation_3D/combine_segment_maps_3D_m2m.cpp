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

#define ARGC_SEGMENT_MAP_1 1
#define ARGC_SEGMENT_MAPPING_1 2
#define ARGC_SEGMENT_MAP_2 3
#define ARGC_SEGMENT_MAPPING_2 4
#define ARGC_N_OVERLAP_SECTION 5
#define ARGC_AREA_OVERLAP_THRESHOLD 6
#define ARGC_AREA_OVERLAP_NORM1_THRESHOLD 6
#define ARGC_AREA_OVERLAP_NORM2_THRESHOLD 6
#define ARGC_OUTPUT_MAPPING_1 9
#define ARGC_OUTPUT_MAPPING_2 10
int main(int argc, char * argv[]){
	if(argc<ARGC_OUTPUT_MAPPING_2){
		printf("Usage combine_segment_maps_3D_aov0_m2m_b segment_stack_1 segment_mapping_1 segment_stack_2 segment_mapping_2 area_overlap_threshold area_overlap_norm1_threshold area_overlap_norm2_threshold output_mapping_1 output_mapping_2\n");
		printf("input params:\n");
		printf("\t1. segment label map 1 MxNxP uint32 matrix\n");
		printf("\t2. segment mapping 1: map labels in segment_stack to merged labels Qx1 uint32\n");
		printf("\t3. segment label map 2 MxNxP uint32 matrix\n");
		printf("\t4. segment mapping 2: map labels in segment_stack to merged labels Qx1 uint32\n");
		printf("\t5. Rule1 area of overlap threshold (pixels)\n");
		printf("\t6. Rule1 normalized area of overlap threshold\n");
		printf("\t7. Rule2 normalized area of overlap threshold\n");
		printf("\t8. output label mapping for segment stack 1\n");
		printf("\t9. output label mapping for segment stack 2\n");
		return 1;
	}

  printf("\nSTART: combine_segment_maps_3D_aov0_m2m_b\n");

  //
  // Read in the first segment stack
  //
	int stack_width_1;
  int stack_height_1;
  int stack_depth_1;
	int stack_n_pixel_1;
  int stack_n_pixel_plane_1;
	unsigned int * segment_map_1;
  {
    int size_element, n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAP_1], size_element, n_dimension,
                                             dimensions);
    if(size_element<=0 || read_data==NULL){
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
    segment_map_1 = (unsigned int *) read_data;
    delete [] dimensions;
  }
  printf("stack_width_1=%d, stack_height_1=%d, stack_depth_1=%d, stack_n_pixel_1=%d, stack_n_pixel_plane_1=%d\n",
            stack_width_1, stack_height_1, stack_depth_1, stack_n_pixel_1, stack_n_pixel_plane_1);

  unsigned int max_segment_label_1 = 0;
  {
    int i;
    unsigned int * l = segment_map_1;
    for(i=0; i<stack_n_pixel_1; i++, l++)
      max_segment_label_1 = MAX(max_segment_label_1, *l);
  }
  
  unsigned int * label_mapping_1;
  unsigned int max_label_mapping_1=0;
  {
    int size_element, n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAPPING_1], size_element, n_dimension,
                                             dimensions);
    if(size_element<=0 || read_data==NULL){
      printf("\nERROR: Could not open label mapping\n");
      return -1;
    }
    if(n_dimension!=1){
      printf("\nERROR: Not a 1D label mapping\n");
      return -1;
    }
    if(dimensions[0]!=max_segment_label_1+1){
      printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
      return -1;
    }
    label_mapping_1 = (unsigned int *) read_data;
    for(int i=0; i<max_segment_label_1; i++)
      max_label_mapping_1 = MAX(max_label_mapping_1, label_mapping_1[i]);
    delete [] dimensions;
  }

  //
  // Read in the second segment stack
  //
	int stack_width_2;
  int stack_height_2;
  int stack_depth_2;
	int stack_n_pixel_2;
  int stack_n_pixel_plane_2;
	unsigned int * segment_map_2;
  {
    int size_element, n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAP_2], size_element, n_dimension,
                                             dimensions);
    if(size_element<=0 || read_data==NULL){
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
    segment_map_2 = (unsigned int *) read_data;
    delete [] dimensions;
  }
  printf("stack_width_2=%d, stack_height_2=%d, stack_depth_2=%d, stack_n_pixel_2=%d, stack_n_pixel_plane_2=%d\n",
            stack_width_2, stack_height_2, stack_depth_2, stack_n_pixel_2, stack_n_pixel_plane_2);
  if(stack_width_1!=stack_width_2 || stack_height_1!=stack_height_2){
    printf("\nERROR: Stack planes do not match.\n");
    return -1;
  }

  unsigned int max_segment_label_2 = 0;
  {
    int i;
    unsigned int * l = segment_map_2;
    for(i=0; i<stack_n_pixel_2; i++, l++)
      max_segment_label_2 = MAX(max_segment_label_2, *l);
  }
  
  unsigned int * label_mapping_2;
  unsigned int max_label_mapping_2=0;
  {
    int size_element, n_dimension, * dimensions;
    void * read_data = iput::io::fread_raw_array(argv[ARGC_SEGMENT_MAPPING_2], size_element, n_dimension,
                                             dimensions);
    if(size_element<=0 || read_data==NULL){
      printf("\nERROR: Could not open label mapping\n");
      return -1;
    }
    if(n_dimension!=1){
      printf("\nERROR: Not a 1D label mapping\n");
      return -1;
    }
    if(dimensions[0]!=max_segment_label_2+1){
      printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
      return -1;
    }
    label_mapping_2 = (unsigned int *) read_data;
    for(int i=0; i<max_segment_label_2; i++)
      max_label_mapping_2 = MAX(max_label_mapping_2, label_mapping_2[i]);
    delete [] dimensions;
  }
  printf("max_segment_label_1:%d, max_segment_label_2:%d\n", max_segment_label_1,
         max_segment_label_2);
  printf("max_label_mapping_1:%d, max_label_mapping_2:%d\n", max_label_mapping_1,
         max_label_mapping_2);
  
  //////////////////////////////////
  // Get the area of overlap
  //////////////////////////////////
  // Get area of labels and overlap
  int n_section_overlap = atoi(argv[ARGC_N_OVERLAP_SECTION]);
  printf("n_section_overlap: %d\n", n_section_overlap);
  std::map<const unsigned int, int, ltstru> label_area_1, label_area_2;
  std::map<const long unsigned int, int, ltstrlu> label_pair_area_overlap;
  unsigned int max_label_overlap_1=0, max_label_overlap_2=0;
  {
    unsigned int *l1 = segment_map_1 + stack_n_pixel_plane_1*(stack_depth_1-n_section_overlap);
    unsigned int *l2 = segment_map_2;
    unsigned int s1, s2;
    int i;
    long unsigned int l_p, l;
    for(i=0; i<n_section_overlap*stack_n_pixel_plane_1; i++, l1++, l2++){
      s1 = label_mapping_1[*l1];
      s2 = label_mapping_2[*l2];
      max_label_overlap_1 = MAX(s1, max_label_overlap_1);
      max_label_overlap_2 = MAX(s2, max_label_overlap_2);
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

  // Get sets of overlapping labels should be linked together.
  // Method: Build a bipartite graph with nodes as segments from
  // the two stacks. Edges are added if two segments are to be
  // linked. Connected components in this graph would indicate
  // sets of linked segments.
  {
    using namespace boost;
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::vertices_size_type size_type;
    int N = max_label_overlap_1 + max_label_overlap_2;
    Graph G(N);

    std::vector<size_type> rank(num_vertices(G));
    std::vector<Vertex> parent(num_vertices(G));
    typedef size_type* Rank;
    typedef Vertex* Parent;
    disjoint_sets<Rank, Parent>  ds(&rank[0], &parent[0]);
    initialize_incremental_components(G, ds);

    // Add label pairs to graph G if they should be linked into one segment
    {
      int area_overlap_threshold = atoi(argv[ARGC_AREA_OVERLAP_THRESHOLD]);
      double area_overlap_norm1_threshold = atof(argv[ARGC_AREA_OVERLAP_THRESHOLD]);
      double area_overlap_norm2_threshold = atof(argv[ARGC_AREA_OVERLAP_THRESHOLD]);
      std::map<const long unsigned int, int, ltstrlu>::iterator pair_it;
      unsigned long int l_p, l;
      for(pair_it=label_pair_area_overlap.begin(); pair_it!=label_pair_area_overlap.end(); pair_it++){
        long unsigned int l_p = (*pair_it).first;
        int area_overlap = (*pair_it).second;
        unsigned int s1, s2;
        s1 = l_p & 0x00000000ffffffff;
        s2 = l_p >> 32;
        if(s1==0 || s2==0)
          continue;
        int area_1 = label_area_1[s1];
        int area_2 = label_area_2[s2];
        bool rule_1 = area_overlap>area_overlap_threshold && area_overlap/MAX(area_1, area_2)>area_overlap_norm1_threshold;
        bool rule_2 = area_overlap/MIN(area_1, area_2)>area_overlap_norm2_threshold;
        bool rule = rule_1 || rule_2;
        if(rule){
          add_edge(s1, s2+max_label_overlap_1, G);
        }
      }
    }

    // Find connected components
    incremental_components(G, ds);

    // Generate combined mapping by modifying the label_mapping's
    unsigned int new_label_offset = max_label_mapping_1 +
      max_label_mapping_2;
    // label_mapping_1
    {
      int i;
      unsigned int l;
      for(i=0; i<max_segment_label_1; i++){
        l = label_mapping_1[i];
        if(l>max_label_overlap_1)
          continue; // didn't occur in overlapping region - do nothing
        label_mapping_1[i] = ds.find_set(l) + new_label_offset;
      }
    }
    // label_mapping_2
    {
      int i;
      unsigned int l;
      for(i=0; i<max_segment_label_2; i++){
        l = label_mapping_2[i];
        if(l>max_label_overlap_2){
          label_mapping_2[i] += max_label_mapping_1;
          continue; // didn't occur in overlapping region - just offset
        }
        label_mapping_2[i] = ds.find_set(l+max_label_overlap_1) +
          new_label_offset;
      }
    }
  }

  // Write new label mapping 1
  {
    int size_element = sizeof(unsigned int);
    int n_dimension = 1;
    int dimensions[1] = {max_segment_label_1+1};
    int err;
    err = iput::io::fwrite_raw_array(argv[ARGC_OUTPUT_MAPPING_1],
                                     size_element, n_dimension, dimensions,
                                     (void *) label_mapping_1);
    if(err!=0){
      printf("ERROR combine_segment_maps_3D_aov0_m2m_b: could not write label mapping 1 [%d]\n", err);
      delete [] segment_map_1;
      delete [] label_mapping_1;
      delete [] segment_map_2;
      delete [] label_mapping_2;
      return -2;
    }
  }  
  // Write new label mapping 2
  {
    int size_element = sizeof(unsigned int);
    int n_dimension = 1;
    int dimensions[1] = {max_segment_label_2+1};
    int err;
    err = iput::io::fwrite_raw_array(argv[ARGC_OUTPUT_MAPPING_2],
                                     size_element, n_dimension, dimensions,
                                     (void *) label_mapping_2);
    if(err!=0){
      printf("ERROR combine_segment_maps_3D_aov0_m2m_b: could not write label mapping 2 [%d]\n", err);
      delete [] segment_map_1;
      delete [] label_mapping_1;
      delete [] segment_map_2;
      delete [] label_mapping_2;
      return -2;
    }
  }  

  delete [] segment_map_1;
  delete [] label_mapping_1;
  delete [] segment_map_2;
  delete [] label_mapping_2;
  printf("\nSTOP: combine_segment_maps_3D_aov0_m2m_b\n");
  return 0;
}
