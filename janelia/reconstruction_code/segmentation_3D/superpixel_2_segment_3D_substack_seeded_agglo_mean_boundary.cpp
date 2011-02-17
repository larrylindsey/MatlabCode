// Superpixel to segment 3D seeded agglomerative mean
// boundary. Combine substack segment.
//
// The seeds correspond to segments along the outer shell of the
// stack. All segments are forced to be members of exactly one of the
// seed-segments.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	06162009	init. code
// v1 06242009  seeded version
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <image_lib.h>

#include <file_input_output.h>
#include <segment_3D_utilities.h>
#include <stack_3D_utilities.h>
#include <merge_sets_h.h>

namespace std{
  using namespace __gnu_cxx;
}

//#define __DEBUG__

extern int stack_width, stack_height, stack_depth, stack_n_pixel, stack_n_pixel_plane;
unsigned int * n_pixel_of_label;
void merge_area(unsigned int s1, unsigned int s2){
  n_pixel_of_label[s1] += n_pixel_of_label[s2];
}

// Hash function from pair of labels to unsigned int
typedef unsigned long int Label_Pair;
struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

struct equli{ bool operator()(const Label_Pair s1, const Label_Pair s2) const
    {return s1==s2;}};
typedef  std::hash_map<Label_Pair, unsigned int, std::hash<Label_Pair>, equli> Label_Pair_Hash;

Label_Pair label_pair_2_id(unsigned int s1, unsigned int s2){
  if(s1>s2){
    Label_Pair lp1;
    lp1 = s2;
    lp1 <<= 32;
    Label_Pair lp;
    lp = s1;
    lp |= lp1;
    return lp;
  }
  Label_Pair lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair lp;
  lp = s2;
  lp |= lp1;
  return lp;
}
std::pair<unsigned int, unsigned int> id_2_label_pair(Label_Pair lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

// Priority queue for merging
struct ltdb{ bool operator()(const double f1, const double f2){ return f1<f2; }};
typedef std::map<const double, Label_Pair, ltdb> Label_Pair_Q;

using namespace boost;

// type of the graph: unsorted sequence of edges, unsorted set of vertices,
// undirected graph, no vertex property and edge property Boundary_Stat
typedef adjacency_list<vecS, vecS, undirectedS, no_property> Graph;

typedef struct{
  unsigned int n_boundary_pixel, sum_boundary_value;
  Label_Pair_Q::iterator merge_q_pos; 
} Boundary_Stat_0;
typedef  std::hash_map<Label_Pair, Boundary_Stat_0, std::hash<Label_Pair>, equli> Label_Pair_Edge_Prop_Hash;


#define ARGC_BOUNDARY_NAME_FORMAT 1
#define ARGC_SUPERPIXEL_NAME_FORMAT 2
#define ARGC_SEGMENT_NAME_FORMAT 3
#define ARGC_OUTPUT_NAME_FORMAT 4
#define ARGC_N_SUBSTACK 5
#define ARGC_SUBSTACK_LIMIT_SEQ 6
#define ARGC_F_THRESHOLD_SEQ (ARGC_SUBSTACK_LIMIT_SEQ+2*atoi(argv[ARGC_N_SUBSTACK]))

int main(int argc, char * argv[]){
  if(argc==0){
    printf("Usage stitch_substack_segment_overlap_aov0_b superpixel_name_format segment_name_format area_overlap_threshold area_overlap_norm1_threshold area_overlap_norm2_threshold output_name_format section_sequence\n");
    printf("input params:\n");
    printf("\t1. Format of boundary stack file names.\n");
    printf("\t2. Format of superpixel file names.\n");
    printf("\t3. Format of segment file names.\n");
    printf("\t4. Format of output mapping file names.\n");
    printf("\t5. Sequence of thresholds.\n");
    return 1;
  }

  printf("\nSTART: superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary_b\n");

  unsigned int n_f_threshold = argc - ARGC_F_THRESHOLD_SEQ;
  printf("n_f_threshold = %d\n", n_f_threshold);
	unsigned char * f_thresholds = (unsigned char *)
    new unsigned char[n_f_threshold];
  {
    unsigned int i;
    printf("f_thresholds: ");
    for(i=0; i<n_f_threshold; i++){
      f_thresholds[i] = (unsigned char) atoi(argv[ARGC_F_THRESHOLD_SEQ+i]);
      printf("%u ", (unsigned int) f_thresholds[i]);
    }
    printf("\n");
  }
  char * output_mapping_name_format =
    argv[ARGC_OUTPUT_NAME_FORMAT];

  printf("Building a joint Region Adjacency Graph for all the substacks.\n");
  unsigned int n_substack = atoi(argv[ARGC_N_SUBSTACK]);
  char filename[10000];

  unsigned int * max_superpixel_labels =
    (unsigned int*) new unsigned int[n_substack];
  // offset to make labels non-overlapping across substacks
  unsigned int * label_mapping_offsets =
    (unsigned int*) new unsigned int[n_substack+1];
  label_mapping_offsets[0] = 0;
  
  Label_Pair_Edge_Prop_Hash edge_prop;
  Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
    segment_sets(NULL);

  unsigned int * superpixel_label_map_prev = NULL;
  Hash_UInt32_UInt32 * label_mappings =
    (Hash_UInt32_UInt32 *) new Hash_UInt32_UInt32[n_substack];
  int stack_n_pixel_prev = -1;
  // Seed a subset of the segments. All other segments are forced to
  // merge with one of the seeeded segments.  If the segments along
  // the outer shell of the stack are defined as seeds then we are
  // forcing all bodies in the final segmentation to have a connection
  // with the outer shell. For this case, all segments touching the
  // outer shell of the stack are labelled to be seeeds.
  Hash_UInt32_UInt32 is_seeded;
  
#ifdef __PROFILE__
#endif
  for(unsigned int substack_id=0; substack_id<n_substack; substack_id++){
    printf("----------------------------\n");
    printf("substack_id: %d\n", substack_id);
    unsigned int substack_start, substack_end;
    substack_start = atoi(argv[ARGC_SUBSTACK_LIMIT_SEQ+substack_id*2]);
    substack_end = atoi(argv[ARGC_SUBSTACK_LIMIT_SEQ+substack_id*2+1]);
    printf("substack_start: %d, substack_end: %d\n", substack_start,
           substack_end);
    
    printf("Reading the boundary map stack.\n");
    Stack * stack;
    {
      sprintf(filename, argv[ARGC_BOUNDARY_NAME_FORMAT], substack_start,
              substack_end);
      printf("filename: %s\n", filename);
      stack = Read_Stack(filename);
      stack_width = stack->width;
      stack_height = stack->height;
      stack_depth = stack->depth;
      stack_n_pixel = stack_height*stack_width*stack_depth;
      stack_n_pixel_plane = stack_width*stack_height;
      printf("stack_width=%d, stack_height=%d, stack_depth=%d, stack_n_pixel=%d, stack_n_pixel_plane=%d\n",
             stack_width, stack_height, stack_depth,
             stack_n_pixel, stack_n_pixel_plane);
    }
    
    printf("Reading in the superpixel map\n");
    unsigned int * superpixel_label_map;
    {
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      sprintf(filename, argv[ARGC_SUPERPIXEL_NAME_FORMAT], substack_start,
              substack_end);
      printf("filename: %s\n", filename);
      void * read_data = iput::io::fread_raw_array(filename,
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
      if(dimensions[0]!=stack_width || dimensions[1]!=stack_height ||
         dimensions[2]!=stack_depth){
        printf("Segment map's dimensions don't match with those of image stack\n");
        printf("Segment map's dimensions: width:%d height:%d depth:%d\n",
               dimensions[0], dimensions[1], dimensions[2]);
        return -1;
      }
      superpixel_label_map = (unsigned int *) read_data;
      for(unsigned int i=0; i<stack_n_pixel; i++)
        max_superpixel_labels[substack_id] =
          std::max(max_superpixel_labels[substack_id],
                   superpixel_label_map[i]);
    }
 
    printf("Reading in the segment label mapping\n");
    {
      iput::io::Type_Element type_element;
      int n_dimension, * dimensions;
      sprintf(filename, argv[ARGC_SEGMENT_NAME_FORMAT], substack_start,
              substack_end);
      printf("filename: %s\n", filename);
      void * read_data = iput::io::fread_raw_array(filename,
                                                   type_element, n_dimension,
                                                   dimensions);
      if(iput::io::is_error_type_element(type_element) || read_data==NULL){
        printf("\nERROR: Could not open label mapping\n");
        return -1;
      }
      if(!(n_dimension==1 || (n_dimension==2 && dimensions[1]==1))){
        printf("\nERROR: Not a 1D label mapping\n");
        return -1;
      }
      if(dimensions[0]!=max_superpixel_labels[substack_id]+1){
        printf("\nERROR: Label mapping's dimensions don't match with max segment label\n");
        return -1;
      }
      {
        unsigned int i;
        unsigned int * label_mapping = (unsigned int *) read_data;
        for(i=0; i<max_superpixel_labels[substack_id]+1; i++)
          label_mappings[substack_id][i] = label_mapping[i];
      }
      
      {
        unsigned int max_label_mapping=0;
        for(int i=0; i<max_superpixel_labels[substack_id]+1; i++){
          label_mappings[substack_id][i] +=
            label_mapping_offsets[substack_id];
          max_label_mapping = std::max(max_label_mapping,
                                       label_mappings[substack_id][i]);
        }
        label_mapping_offsets[substack_id+1] = max_label_mapping + 1;
      }
      
      delete [] dimensions;
    }
  
    unsigned char * boundary_value = (unsigned char *) stack->array;
 
    printf("Setting segments along outer shell to be seeds ...\n");
    {
      unsigned int i,j;
      if(substack_id==0){
        // top face. Skip the bottom face as it is overlapping
        for(i=0; i<stack_n_pixel_plane; i++)
          is_seeded[label_mappings[substack_id][superpixel_label_map[i]]] = 1;
      }
      if(substack_id==n_substack-1){
        // bottom face.
        for(i=stack_n_pixel-stack_n_pixel_plane; i<stack_n_pixel; i++)
          is_seeded[label_mappings[substack_id][superpixel_label_map[i]]] = 1;
      }
      for(i=0; i<stack_depth; i++){
        for(j=i*stack_n_pixel_plane; j<i*stack_n_pixel_plane+stack_width;
            j++){
          // north face
          is_seeded[label_mappings[substack_id][superpixel_label_map[j]]] = 1;
          // south face
          is_seeded[label_mappings[substack_id][
              superpixel_label_map[j + stack_width*(stack_height-1)]]] = 1;
        }
        for(j=i*stack_n_pixel_plane;
            j<i*stack_n_pixel_plane+stack_width*(stack_height-1);
            j+=stack_width){
          // west face
          is_seeded[label_mappings[substack_id][superpixel_label_map[j]]] = 1;
          // south face
          is_seeded[label_mappings[substack_id][
              superpixel_label_map[j + stack_width - 1]]] = 1;
        }
      }
    }
    
    printf("Collecting label-pair statistics within substack ...\n");
    {
      // skip the overlapping planes at the end
      int stack_depth_orig = stack_depth;
      if(substack_id!=n_substack-1)
        stack_depth -= (int)substack_end -
          atoi(argv[ARGC_SUBSTACK_LIMIT_SEQ+(1+substack_id)*2]) + 1;
      stack_n_pixel = stack_n_pixel_plane * stack_depth;
      
      unsigned int i;
      unsigned int neighbors[27];
      int n_neighbor, n1, n2;
      unsigned int * s = superpixel_label_map;
      unsigned int s1, s2;
      unsigned char * b = boundary_value;
      bool is_boundary;
      for(i=0; i<stack_n_pixel; i++, s++, b++){
        if(i%10000000==0)
          printf("i: %d\n", i);
        s1 = label_mappings[substack_id][*s];
        if(s1!=0){
          n_neighbor = get_neighbor_indexes_6(i, neighbors);
          for(n2=0; n2<n_neighbor; n2++){
            s2 = label_mappings[substack_id][
              superpixel_label_map[neighbors[n2]]];
            if(s2==0)
              continue;
            if(s1==s2)
              continue;
            Label_Pair lp = label_pair_2_id(s1, s2);
            edge_prop[lp].sum_boundary_value += *b;
            edge_prop[lp].n_boundary_pixel ++;
          }
        }

        segment_sets.add_new_set_inc(s1);
        
#ifdef __DEBUG__
          printf("segment_sets.adam[%u]: %u\n", s1, segment_sets.get_adam(s1));
#endif //__DEBUG__
      }

      stack_depth = stack_depth_orig;
    }
    
    if(substack_id!=0){
      printf("Collecting label-pair statistics between this and\n");
      printf("the previous substack ...\n");
      unsigned char * b = boundary_value;
      unsigned int i;
      unsigned int offset_prev = stack_n_pixel_prev - stack_n_pixel_plane;
      for(i=0; i<stack_n_pixel_plane; i++, b++){
        unsigned int s1 = label_mappings[substack_id-1][
          superpixel_label_map_prev[i+offset_prev]];
        unsigned int s2 = label_mappings[substack_id][
          superpixel_label_map[i]];
        Label_Pair lp = label_pair_2_id(s1, s2);
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
      }
    }

    Kill_Stack(stack);
    if(substack_id!=0){
      delete [] superpixel_label_map_prev;
    }
    superpixel_label_map_prev = superpixel_label_map;
    stack_n_pixel_prev = stack_n_pixel;
    
    printf("----------------------------\n");
#ifdef __PROFILE__
#endif
  }
  
  if(superpixel_label_map_prev!=NULL){
    delete [] superpixel_label_map_prev;
  }

  printf("Building Region Adjacency Graph ...\n");
  // Construct region adjacency graph with segments as nodes and edges
  // between spatially adjacent segments.  create a typedef for the
  // Graph type
#ifdef __PROFILE__
#endif
  Graph g(10000);
  Label_Pair_Hash label_adjacency;
  Label_Pair_Q merge_q;
  {
    Label_Pair_Q::iterator q_it;
    Label_Pair_Edge_Prop_Hash::iterator lp_it;
    unsigned int s0, s1;
    unsigned int n, b;
    double m;
    Label_Pair lp;
    unsigned int i=0;
    for(lp_it=edge_prop.begin(); lp_it!=edge_prop.end();
        lp_it++){
      if(i%10000==0)
        printf("i: %d\n", i);
      i++;
      lp = (*lp_it).first;
      tie(s0,s1) = id_2_label_pair(lp);
      if(s0==0 || s1==0 || s0==s1)
        continue;
      graph_traits<Graph>::edge_descriptor e;
      bool is_inserted;
      label_adjacency[lp] = 1;
      tie(e, is_inserted) = add_edge(s0, s1, g);
      n = (*lp_it).second.n_boundary_pixel;
      b = (*lp_it).second.sum_boundary_value;
#ifdef __DEBUG__
      printf("s0: %u, s1: %u, sum_boundary_value: %u, n_boundary_pixel: %u\n",
             s0, s1, b, n);
#endif //_DEBUG__
      m = (double)b/(double)n; // mean boundary value
      // insert in merge queue if exactly one of the segments is seeded.
      if((is_seeded[s0]>0) ^ (is_seeded[s1]>0)){
        q_it = merge_q.find(m); // keep priority keys unique.
        while(q_it!=merge_q.end()){
          m += 0.00001*((double)(rand()%1000)); // epsilon increments
          q_it = merge_q.find(m);
        }
        {
          const std::pair<double, Label_Pair> v(m, lp);
          bool is_inserted;
          tie(q_it, is_inserted) = merge_q.insert(v);
        }
        (*lp_it).second.merge_q_pos = q_it;
      }
    }
  }
#ifdef __PROFILE__
#endif

  printf("Performing agglo. segmentation ...\n");
  // Make mergers. Save label mapping for f_threshold_seq
  {
    Label_Pair_Q::iterator q_it, q_it1;
    double m, m1; // merge criterion for a pair
    unsigned int s0, s1; // segments to be merged
    unsigned int n0, n1; // their respective neighbors
    int f_threshold_id=0;
    double f_threshold;
    graph_traits<Graph>::out_edge_iterator edge_0, edge_0_end,
      edge_1, edge_1_end; // to iterate over the neighbors of s0, s1
    graph_traits<Graph>::edge_descriptor e0, e1;
    unsigned int b, n; // boundary value and number of pixels

    // get number of neighbors for each vertex
    Hash_UInt32_UInt32 n_neighbor;
    {
      graph_traits<Graph>::vertex_iterator vs, ve;
      tie(vs, ve) = vertices(g);
      for(; vs!=ve; vs++)
        n_neighbor[(*vs)] = out_degree(*vs, g);
    }
    q_it = merge_q.begin();
    if(q_it==merge_q.end())
      return 0;
    m = (*q_it).first;
    tie(s0,s1) = id_2_label_pair((*q_it).second);
    while(q_it!=merge_q.end() && f_threshold_id<n_f_threshold){
      f_threshold = f_thresholds[f_threshold_id];
      printf("f_threshold:%g, m:%g\n", f_threshold, m);
#ifdef __PROFILE__
#endif
      while(q_it!=merge_q.end() && m<f_threshold){ // merge segments s0 and s1
        if(s0==0 || s1==0 || (is_seeded[s0]>0 && is_seeded[s1]>0)){
          merge_q.erase(q_it);
          q_it = merge_q.begin();
          if(q_it==merge_q.end())
            break;
          m = (*q_it).first;
          tie(s0,s1) = id_2_label_pair((*q_it).second);
          
          f_threshold = f_thresholds[f_threshold_id];
          continue;
        }
        
        printf("m: %g, s0: %u, s1:%u\n", m, s0, s1);
        // merge the superpixel sets of s1 into s0
        segment_sets.merge(s0,s1);

        // remove edge from graph
        label_adjacency[label_pair_2_id(s0,s1)]=0;
        
        // remove from merge queue
        merge_q.erase(q_it);

        // merge the neighbors of segment s1 into s0. Exactly one of s0
        // and s1 is seeded. If s1 is seeded then swap them.
        if(is_seeded[s1]>0){
          unsigned int s01 = s0;
          s0 = s1;
          s1 = s01;
        }
        // From now s0 is seeded and s1 is not seeded.
        
        n_neighbor[s0] += n_neighbor[s1];
        
        tie(edge_1, edge_1_end) = out_edges(s1, g);
        while(edge_1!=edge_1_end){
          n1 = target(*edge_1, g);
          Label_Pair lp_n1_s1 = label_pair_2_id(n1,s1);
          if(label_adjacency[lp_n1_s1]==0){
            edge_1++;
            continue;
          }
          
          Label_Pair lp_n1_s0 = label_pair_2_id(n1,s0);
          if(n_neighbor[n1]!=0){
            if(label_adjacency[lp_n1_s0]==1){

              label_adjacency[lp_n1_s1] = 0;

              Label_Pair_Edge_Prop_Hash::iterator ep_it_n1_s0 =
                edge_prop.find(lp_n1_s0);
              Label_Pair_Edge_Prop_Hash::iterator ep_it_n1_s1 =
                edge_prop.find(lp_n1_s1);
              
              n = (*ep_it_n1_s0).second.n_boundary_pixel +=
                (*ep_it_n1_s1).second.n_boundary_pixel;
              b = (*ep_it_n1_s0).second.sum_boundary_value +=
                (*ep_it_n1_s1).second.sum_boundary_value;
              
              // erase previous entry in merge Q.
              if(is_seeded[n1]>0) // exactly one of n1 and s1 was seeded
                (*((*ep_it_n1_s1).second.merge_q_pos)).second = 0;
              else // exactly one of n1 and s0 was seeded
                (*((*ep_it_n1_s0).second.merge_q_pos)).second = 0;

              if(is_seeded[n1]==0){ // n1 is not seeded so add to merge Q
                // insert in merge queue
                // random increment for more likely uniqueness
                m1 = (double)b/(double)n  +
                  0.000001*((double)(rand()%9973));
                q_it1 = merge_q.find(m1);
                while(q_it1!=merge_q.end()){
                  m1 += 0.000001*((double)(rand()%9973)); // epsilon increments
                  q_it1 = merge_q.find(m1);
                }
                {
                  const std::pair<double, Label_Pair>
                    v(m1, lp_n1_s0);
                  bool is_inserted;
                  tie(q_it1, is_inserted) = merge_q.insert(v);
                }
                // save the new entry's Q pos
                (*ep_it_n1_s0).second.merge_q_pos = q_it1;
              }
            }
            else{
              // move this neighbor of s1 to s0
              graph_traits<Graph>::edge_descriptor e;
              bool is_inserted;
              tie(e, is_inserted) = add_edge(s0, n1, g);
              
              Label_Pair_Edge_Prop_Hash::iterator ep_it_n1_s1 =
                edge_prop.find(lp_n1_s1);
              
              label_adjacency[lp_n1_s0] = 1;

              Boundary_Stat_0 e_p;
              e_p.n_boundary_pixel = (*ep_it_n1_s1).second.n_boundary_pixel;
              e_p.sum_boundary_value = (*ep_it_n1_s1).second.sum_boundary_value;
              // Regarding merge Q for this neighbor
              if(is_seeded[n1]>0){
                // n1 was seeded. Since s1 is not seeded, remove the Q
                // node <n1,s1>
                (*((*ep_it_n1_s1).second.merge_q_pos)).second = 0;
              }
              else{
                // There is no Q node <n1,s1>. Therefore, add <n1,s0> 
                // insert in merge queue
                // random increment for more likely uniqueness
                m1 = (double)e_p.sum_boundary_value /
                  (double)e_p.n_boundary_pixel  +
                  0.000001*((double)(rand()%9973));
                q_it1 = merge_q.find(m1);
                while(q_it1!=merge_q.end()){
                  // epsilon increments
                  m1 += 0.000001*((double)(rand()%9973));
                  q_it1 = merge_q.find(m1);
                }
                {
                  const std::pair<double, Label_Pair>
                    v(m1, lp_n1_s0);
                  bool is_inserted;
                  tie(q_it1, is_inserted) = merge_q.insert(v);
                }
                // save the new entry's Q pos
                e_p.merge_q_pos = q_it1;
              }
              
              edge_prop[lp_n1_s0] = e_p;
            }
          }

          label_adjacency[lp_n1_s1]=0;
          
          edge_1++;
        }

        // remove s1
        //clear_vertex(s1, g);
        n_neighbor[s1] = 0;
        
        q_it = merge_q.begin();
        if(q_it==merge_q.end())
          break;
        m = (*q_it).first;
        tie(s0,s1) = id_2_label_pair((*q_it).second);
        
        f_threshold = f_thresholds[f_threshold_id];
        //printf("f_threshold:%g, m:%g\n", f_threshold, m);
      }
      
      // copy the current superpixel-to-segment mapping into
      // output variable and save to file
      printf("Dumping new label mapping ...\n");
      for(unsigned int substack_id=0; substack_id<n_substack; substack_id++){
        printf("substack_id: %u\n", substack_id);
        unsigned int substack_start =
          atoi(argv[ARGC_SUBSTACK_LIMIT_SEQ+substack_id*2]);
        unsigned int substack_end =
          atoi(argv[ARGC_SUBSTACK_LIMIT_SEQ+substack_id*2+1]);
        printf("substack_start: %u, substack_end: %u\n",
               substack_start, substack_end);
        int err;
        int size_element = sizeof(unsigned int);
        iput::io::Type_Element type_element = iput::io::get_type_element(
          iput::io::type_id_unsigned_int(), size_element);
        int n_dimension = 1;
        int dimensions[1] = {max_superpixel_labels[substack_id]+1};

        unsigned int  * sp_2_seg_map = (unsigned int *) new
          unsigned int[max_superpixel_labels[substack_id]+1];
        for(unsigned int i=0; i<max_superpixel_labels[substack_id]+1; i++)
          sp_2_seg_map[i]=0;
        
        Hash_UInt32_UInt32::iterator l_it;
#ifdef __DEBUG__
        printf("New label mapping:\n");
#endif // __DEBUG__
        for(l_it=label_mappings[substack_id].begin();
            l_it!=label_mappings[substack_id].end(); l_it++){
          sp_2_seg_map[(*l_it).first] = segment_sets.get_adam((*l_it).second);
#ifdef __DEBUG__
          printf("segment_sets.adam[%u]: %u\n", (*l_it).first,
                 segment_sets.get_adam((*l_it).second));
#endif // __DEBUG__
        }
        char * file_name = (char *) new char[strlen(output_mapping_name_format)
                                             + 30];
        sprintf(file_name, output_mapping_name_format,
                substack_start, substack_end, (unsigned int) f_threshold);
        err = iput::io::fwrite_raw_array(file_name, type_element,
                                         n_dimension, dimensions,
                                         (void *) sp_2_seg_map);
        delete [] file_name;
        delete [] sp_2_seg_map;

        if(err!=0){
          printf("ERROR superpixel_2_segment_3D_ladder_b: could not write label mapping [%d]\n", err);
          return -2;
        }
      }        
#ifdef __PROFILE__
#endif
      f_threshold_id++;
    }
  }
  
  printf("\nSTOP: superpixel_2_segment_3D_substack_seeded_agglo_mean_boundary_b\n");
  return 0;
}
