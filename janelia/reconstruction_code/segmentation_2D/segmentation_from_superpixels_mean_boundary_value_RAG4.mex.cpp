// compute segmentation from a superpixel map (agglomerative
// clustering) Feature: between each segment, compute the mean
// boundary value and merge them if this value is below a specified
// threshold
// 
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	08042008	init. code
//

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <merge_sets.h>
#include <merge_sets_h.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

namespace std{
  using namespace __gnu_cxx;
}

struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
  {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

// Hash function from pair of labels to unsigned int
typedef unsigned long int Label_Pair;
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

typedef std::hash_map<Label_Pair, Label_Pair_Q::iterator, std::hash<Label_Pair>, equli> Hash_Label_Pair_Q_It;

typedef struct{
  unsigned int n_boundary_pixel, sum_boundary_value;
  Label_Pair_Q::iterator merge_q_pos; 
} Boundary_Stat_0;
typedef  std::hash_map<Label_Pair, Boundary_Stat_0, std::hash<Label_Pair>, equli> Label_Pair_Edge_Prop_Hash;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage old_label_to_new_label = segmentation_from_superpixels_mean_boundary_value_RAG(superpixel_label_map, boundary, f_threshold_seq, boundary_length_seq);\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. superpixel label map MxN matrix (uint32)\n");
    mexPrintf("\t2. scalar field on which watershed was performed MxN matrix (uint8)\n");
    mexPrintf("\t3. f_thresholds for which segmentation is to be saved 1xR (double).\n");
    mexPrintf("\t4. minimum length of boundaries to be considered for merging 1xR (uint32).\n");    
    mexPrintf("output:\n");
    mexPrintf("\t1. mapping from input segment labels to merged segment labels 1xR cell array of 1xS uint32 (S=max. label id).\n");
    return;
  }
  if(nrhs!=4){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: segmentation_from_superpixels_mean_boundary_value_RAG\n");
  int numDim = mxGetNumberOfDimensions(prhs[0]);
  if(numDim!=2){
    mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
    return;
  }

  const int * sizeImage;
  sizeImage = mxGetDimensions(prhs[0]);
  int width = sizeImage[1], height = sizeImage[0];
  int n_pixel = height*width;
  unsigned int * superpixel_label_map = (unsigned int *) mxGetPr(prhs[0]);
	
  unsigned char * boundary_map = (unsigned char *) mxGetPr(prhs[1]);
  double * f_thresholds = mxGetPr(prhs[2]);
  int n_f_threshold = mxGetN(prhs[2]);
  unsigned int * boundary_length_thresholds = (unsigned int *) mxGetPr(prhs[3]);

  // get sum and number of boundary values for each pair of adjacent segments
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
  unsigned int max_superpixel_id = 0;
  Label_Pair_Hash label_adjacency;
  Label_Pair_Edge_Prop_Hash edge_prop;
  Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
    superpixel_sets(NULL);
  {
    int x, y;
    unsigned int * s;
    unsigned char * b;
    Label_Pair lp;
    Label_Pair_Edge_Prop_Hash::iterator ep_it;
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
      }
      else{
        max_superpixel_id = MAX(max_superpixel_id, *s);
        superpixel_sets.add_new_set_inc(*s);
      }
    }
    for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
      }
      else{
        max_superpixel_id = MAX(max_superpixel_id, *s);
        superpixel_sets.add_new_set_inc(*s);
      }
    }
    for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
      }
      else{
        max_superpixel_id = MAX(max_superpixel_id, *s);
        superpixel_sets.add_new_set_inc(*s);
      }
    }
    for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
      }
      else{
        max_superpixel_id = MAX(max_superpixel_id, *s);
        superpixel_sets.add_new_set_inc(*s);
      }
    }
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height; // skip first column
        x<width-1;
        x++){
      s++; // skip first element in column
      b++;
      for(y=1;
          y<height-1;
          y++, s++, b++){
        if(*s==0){
          lp = label_pair_2_id(*(s-height), *(s+height));
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
        }
        else{
          max_superpixel_id = MAX(max_superpixel_id, *s);
          superpixel_sets.add_new_set_inc(*s);
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-1), *(s+1));
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height-1), *(s+height+1));
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height+1), *(s+height-1));
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
        }
      }
      s++; // skip last element in column.
      b++;
    }
#ifdef DEBUG_
    {
      Label_Pair_Edge_Prop_Hash::iterator lp_it;
      unsigned int s0, s1;
      Label_Pair lp;
      for(lp_it=edge_prop.begin(); lp_it!=edge_prop.end(); lp_it++){
        lp = (*lp_it).first;
        tie(s0,s1) = id_2_label_pair(lp);
        mexPrintf(
		  "label0: %u, label1: %u, sum_boundary_value: %u, n_boundary_pixel: %u\n",
		  s0, s1, (*lp_it).second.sum_boundary_value,
		  (*lp_it1).second.n_boundary_pixel);
      }
    }
#endif
  }
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif
    
  // Construct region adjacency graph with segments as nodes and edges between
  // spatially adjacent segments.
  // create a typedef for the Graph type
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
  Graph g(MAX(10000, width*height/80));
  Label_Pair_Q merge_q;
  Hash_Label_Pair_Q_It merge_q_pos;
  {
    Label_Pair_Q::iterator q_it;
    Label_Pair_Edge_Prop_Hash::iterator lp_it;
    unsigned int s0, s1;
    unsigned int n, b;
    double m;
    Label_Pair lp;
    for(lp_it=edge_prop.begin(); lp_it!=edge_prop.end();
        lp_it++){
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
      m = (double)b/(double)n; // mean boundary value
      // insert in merge queue
      q_it = merge_q.find(m); // keep priority keys unique.
      while(q_it!=merge_q.end()){
        m += 0.00001*((double)(rand()%1000)); // epsilon increments
        q_it = merge_q.find(m);
      }
      {
        const std::pair<double, Label_Pair> v(m, lp);
        bool is_inserted;
        tie(q_it, is_inserted)  = merge_q.insert(v);
      }
      (*lp_it).second.merge_q_pos = q_it;
    }
  }
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif

  // Make mergers. Save label mapping for f_threshold_seq
  //  Merge_Sets superpixel_sets1(max_superpixel_id, NULL);
  plhs[0] = mxCreateCellMatrix(1, n_f_threshold);
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
    m = (*q_it).first;
    tie(s0,s1) = id_2_label_pair((*q_it).second);
    while(f_threshold_id<n_f_threshold){
      f_threshold = f_thresholds[f_threshold_id];
      mexPrintf("f_threshold:%g, m:%g\n", f_threshold, m);
#ifdef __PROFILE__
      mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
      while(m<f_threshold){ // merge segments s0 and s1
        if(s0==0 || s1==0){
          merge_q.erase(q_it);
          q_it = merge_q.begin();
          m = (*q_it).first;
          tie(s0,s1) = id_2_label_pair((*q_it).second);
          
          f_threshold = f_thresholds[f_threshold_id];
          continue;
        }
        
        // merge the superpixel sets of s1 into s0
        superpixel_sets.merge(s0,s1);

        // remove edge from graph
        label_adjacency[label_pair_2_id(s0,s1)]=0;
        
        // remove from merge queue
        merge_q.erase(q_it);
        
        // merge the neighbors of segments s0 and s1
        if(n_neighbor[s0]<n_neighbor[s1]){
          unsigned int s01 = s0;
          s0 = s1;
          s1 = s01;
        }
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
              // random increment for more likely uniqueness
              m1 = (double)b/(double)n  +
                0.000001*((double)(rand()%9973));
              
              // erase previous entry in merge Q.
              (*((*ep_it_n1_s1).second.merge_q_pos)).second = 0;
              (*((*ep_it_n1_s0).second.merge_q_pos)).second = 0;
              // insert in merge queue
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
              (*((*ep_it_n1_s1).second.merge_q_pos)).second = lp_n1_s0;
              e_p.merge_q_pos = (*ep_it_n1_s1).second.merge_q_pos;
              edge_prop[lp_n1_s0] = e_p;
            }
          }
          
          edge_1++;
        }
        // remove s1
        //clear_vertex(s1, g);
        n_neighbor[s1] = 0;
        
        q_it = merge_q.begin();
        m = (*q_it).first;
        tie(s0,s1) = id_2_label_pair((*q_it).second);
        
        f_threshold = f_thresholds[f_threshold_id];
        //mexPrintf("f_threshold:%g, m:%g\n", f_threshold, m);
      }
      
      // copy the current superpixel-to-segment mapping into
      // output variable
      //      superpixel_sets1.update_adams();
      mxArray * sp_2_seg_map_mx = mxCreateNumericMatrix(
							1, max_superpixel_id+1, mxUINT32_CLASS, mxREAL);
      {
        int i;
        unsigned int * sp_2_seg_map = (unsigned int*) mxGetPr(sp_2_seg_map_mx);
        for(i=0; i<max_superpixel_id+1; i++)
          sp_2_seg_map[i] = superpixel_sets.get_adam(i);
      }
      mxSetCell(plhs[0], f_threshold_id, sp_2_seg_map_mx);
      
#ifdef __PROFILE__
      mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif
      f_threshold_id++;
    }

  }
  mexPrintf("STOP: segmentation_from_superpixels_mean_boundary_value_RAG\n");
  return;
}
