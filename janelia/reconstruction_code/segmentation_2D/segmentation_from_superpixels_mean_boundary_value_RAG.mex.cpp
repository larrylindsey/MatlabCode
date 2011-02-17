// compute segmentation from a superpixel map (agglomerative clustering)
// Feature: between each segment, compute the mean boundary value and merge them if this value is below
// a specified threshold
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
#include <stdlib.h>
#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <merge_sets.h>

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
    unsigned int s=s2;
    s2=s1;
    s1=s;
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

// Edge properties: n_boundary_pixel, sum_boundary_value, merge_queue_position
typedef struct{
  unsigned int n_boundary_pixel, sum_boundary_value;
  Label_Pair_Q::iterator merge_q_pos; 
} Boundary_Stat_0;
enum edge_boundary_t {edge_boundary};
namespace boost{
  BOOST_INSTALL_PROPERTY(edge, boundary);
}
typedef property<edge_boundary_t, Boundary_Stat_0> Boundary_Stat;
// type of the graph: sorted sequence of edges, unsorted set of vertices,
// undirected graph, no vertex property and edge property Boundary_Stat
typedef adjacency_list<setS, vecS, undirectedS, no_property, Boundary_Stat > Graph;

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
  Label_Pair_Hash sum_boundary_value;
  Label_Pair_Hash n_boundary_pixel;
  {
    int x, y;
    unsigned int * s;
    unsigned char * b;
    Label_Pair lp;
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
      else
        max_superpixel_id = MAX(max_superpixel_id, *s);
    }
    for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
      else
        max_superpixel_id = MAX(max_superpixel_id, *s);
    }
    for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
      else
        max_superpixel_id = MAX(max_superpixel_id, *s);
    }
    for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        sum_boundary_value[lp] += *b;
        n_boundary_pixel[lp] ++;
      }
      else
        max_superpixel_id = MAX(max_superpixel_id, *s);
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
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        else
          max_superpixel_id = MAX(max_superpixel_id, *s);
        if(*s==0){
          lp = label_pair_2_id(*(s-1), *(s+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height-1), *(s+height+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        if(*s==0){
          lp = label_pair_2_id(*(s-height+1), *(s+height-1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
      }
      s++; // skip last element in column.
      b++;
    }
#ifdef DEBUG_
    {
      Label_Pair_Hash::iterator lp_it, lp_it1;
      unsigned int s0, s1;
      Label_Pair lp;
      for(lp_it=sum_boundary_value.begin(), lp_it1=n_boundary_pixel.begin();
          lp_it!=sum_boundary_value.end(); lp_it++, lp_it1++){
        lp = (*lp_it).first;
        tie(s0,s1) = id_2_label_pair(lp);
        mexPrintf(
          "label0: %u, label1: %u, sum_boundary_value: %u, n_boundary_pixel: %u\n",
          s0, s1, (*lp_it).second, (*lp_it1).second);
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
  property_map<Graph, edge_boundary_t>::type boundary_stats = get(edge_boundary_t(), g);
  Label_Pair_Q merge_q;
  {
    Label_Pair_Q::iterator q_it;
    Label_Pair_Hash::iterator lp_it, lp_it1;
    unsigned int s0, s1;
    unsigned int n, b;
    double m;
    Label_Pair lp;
    for(lp_it=sum_boundary_value.begin(); lp_it!=sum_boundary_value.end();
        lp_it++){
      lp = (*lp_it).first;
      tie(s0,s1) = id_2_label_pair(lp);
      if(s0==0 || s1==0 || s0==s1)
        continue;
      graph_traits<Graph>::edge_descriptor e;
      bool is_inserted;
      tie(e, is_inserted) = add_edge(s0, s1, g);
      n = n_boundary_pixel[lp];
      b = (*lp_it).second;
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
      boundary_stats[e].n_boundary_pixel = n;
      boundary_stats[e].sum_boundary_value = b;
      boundary_stats[e].merge_q_pos = q_it;
    }
  }
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif

  // Make mergers. Save label mapping for f_threshold_seq
  Merge_Sets superpixel_sets(max_superpixel_id, NULL);
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
        remove_edge(s0, s1, g);
        // remove from merge queue
        merge_q.erase(q_it);
        
        // merge the neighbors of segments s0 and s1
        if(n_neighbor[s0]<n_neighbor[s1]){
          unsigned int s01 = s0;
          s0 = s1;
          s1 = s01;
        }
        n_neighbor[s0] += n_neighbor[s1];
        
        tie(edge_0, edge_0_end) = out_edges(s0, g);
        tie(edge_1, edge_1_end) = out_edges(s1, g);
        while(edge_1!=edge_1_end){
          if(edge_0!=edge_0_end)
            n0 = target(*edge_0, g);
          else
            n0 = 0;
          n1 = target(*edge_1, g);
          
          if(n0!=0 && n0==n1){
            if(n_neighbor[n0]!=0){
              e0 = *edge_0;
              e1 = *edge_1;
              n = boundary_stats[e0].n_boundary_pixel +=
                boundary_stats[e1].n_boundary_pixel;
              b = boundary_stats[e0].sum_boundary_value +=
                boundary_stats[e1].sum_boundary_value;
              m1 = (double)b/(double)n;
              // erase previous entry in merge Q.
//              merge_q.erase(boundary_stats[e0].merge_q_pos);
//              merge_q.erase(boundary_stats[e1].merge_q_pos);
              (*boundary_stats[e0].merge_q_pos).second = 0;
              (*boundary_stats[e1].merge_q_pos).second = 0;
              // insert in merge queue
              q_it1 = merge_q.find(m1); // keep priority keys unique.
              while(q_it1!=merge_q.end()){
                m1 += 0.00002*((double)(rand()%1000)); // epsilon increments
                q_it1 = merge_q.find(m1);
              }
              {
                const std::pair<double, Label_Pair>
                  v(m1, label_pair_2_id(s0,n0));
                bool is_inserted;
                tie(q_it1, is_inserted) = merge_q.insert(v);
              }
              // save the new entry's Q pos
              boundary_stats[e0].merge_q_pos = q_it1;
            }
            edge_0++;
            edge_1++;
          }
          else{
            if(n0>n1 || edge_0==edge_0_end){
              if(n_neighbor[n1]!=0){
                // move this neighbor of s1 to s0
                e1 = graph_traits<Graph>::edge_descriptor(*edge_1);            
                graph_traits<Graph>::edge_descriptor e;
                bool is_inserted;
                tie(e, is_inserted) = add_edge(s0, n1, g);
                boundary_stats[e].n_boundary_pixel =
                  boundary_stats[e1].n_boundary_pixel;
                boundary_stats[e].sum_boundary_value =
                  boundary_stats[e1].sum_boundary_value;
                /*
                n = boundary_stats[e].n_boundary_pixel;
                b = boundary_stats[e].sum_boundary_value;
                m1 = (double)b/(double)n;
                
                // erase previous entry in merge Q.
                merge_q.erase(boundary_stats[e1].merge_q_pos);
                // insert in merge queue
                q_it1 = merge_q.find(m1); // keep priority keys unique.
                while(q_it1!=merge_q.end()){
                m1 += 0.000001; // epsilon increments
                q_it1 = merge_q.find(m1);
                }
                {
                const std::pair<double, Label_Pair>
                v(m1, label_pair_2_id(s0, n1));
                bool is_inserted;
                tie(q_it1, is_inserted) = merge_q.insert(v);
                }
                boundary_stats[e].merge_q_pos = q_it1;
                */
                (*boundary_stats[e1].merge_q_pos).second =
                  label_pair_2_id(s0, n1);
                boundary_stats[e].merge_q_pos = boundary_stats[e1].merge_q_pos;
              }
              edge_1++;
            }
            else{
              if(edge_0!=edge_0_end && n0<n1) // do nothing
                edge_0++;
            }
          }
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
      superpixel_sets.update_adams();
      mxArray * sp_2_seg_map_mx = mxCreateNumericMatrix(
        1, max_superpixel_id+1, mxUINT32_CLASS, mxREAL);
      {
        int i;
        unsigned int * sp_2_seg_map = (unsigned int*) mxGetPr(sp_2_seg_map_mx);
        for(i=0; i<max_superpixel_id+1; i++)
          sp_2_seg_map[i] = superpixel_sets.adam[i];
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
