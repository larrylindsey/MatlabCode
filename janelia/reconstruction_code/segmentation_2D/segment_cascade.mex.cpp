// Apply a cascade of classifiers for agglomerative segmentation
// 
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	01212010  init code
//

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ext/hash_map>
#include <string>
#include <sstream>
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
typedef std::pair<unsigned int,unsigned int> Label_Pair;
typedef unsigned long int Label_Pair_Id;
struct equli{ bool operator()(const Label_Pair_Id s1, const Label_Pair_Id s2) const
  {return s1==s2;}};
typedef  std::hash_map<Label_Pair_Id, unsigned int, std::hash<Label_Pair_Id>, equli> Label_Pair_Id_Hash;

Label_Pair_Id label_pair_2_id(unsigned int s1, unsigned int s2){
  if(s1>s2){
    Label_Pair_Id lp1;
    lp1 = s2;
    lp1 <<= 32;
    Label_Pair_Id lp;
    lp = s1;
    lp |= lp1;
    return lp;
  }
  Label_Pair_Id lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair_Id lp;
  lp = s2;
  lp |= lp1;
  return lp;
}
std::pair<unsigned int, unsigned int> id_2_label_pair(Label_Pair_Id lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

#define N_DIM_BOUNDARY_FEATURE 32
typedef struct{
  unsigned int n_boundary_pixel, sum_boundary_value;
  int boundary_feature[N_DIM_BOUNDARY_FEATURE];
} Boundary_Stat_0;
typedef  std::hash_map<Label_Pair_Id, Boundary_Stat_0, std::hash<Label_Pair_Id>, equli> Label_Pair_Id_Edge_Prop_Hash;

typedef unsigned char BOUNDARY_VALUE;

void num_boundary_pixels_b_less_than(char * stage_type, char * stage_param,
                                     unsigned int * superpixel_label_map,
                                     BOUNDARY_VALUE * boundary_map,
                                     int width, int height,
                                     Label_Pair_Id_Edge_Prop_Hash & edge_prop,
				     Label_Pair_Id_Hash & to_merge_pairs,
				     bool initialize_to_merge_pairs);
void num_boundary_pixels_b_greater_than(char * stage_type, char * stage_param,
					unsigned int * superpixel_label_map,
					BOUNDARY_VALUE * boundary_map,
					int width, int height,
				Label_Pair_Id_Edge_Prop_Hash & edge_prop,
					Label_Pair_Id_Hash & to_merge_pairs,
					bool initialize_to_merge_pairs);

void perform_mergers(unsigned int * label_map, int n_pixel, Label_Pair_Id_Hash & to_merge_pairs);

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    if(nlhs==1){
      plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
      *mxGetPr(plhs[0]) = 1;
      return;
    }
    mexPrintf("Usage old_label_to_new_label = segment_cascade(initial_segment_map, boundary_map, segment_cascade_config\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. starting segment label map MxN matrix (uint32)\n");
    mexPrintf("\t2. scalar field of boundary confidences MxN matrix (uint8)\n");
    mexPrintf("\t3. cascade configuration Rx2 cell array\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. mapping from input segment labels to merged segment labels Sx2 matrix (uint32).\n");
    return;
  }
  if(nrhs!=3){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: segment_cascade\n");
  int numDim = mxGetNumberOfDimensions(prhs[0]);
  if(numDim!=2){
    mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
    return;
  }

  const mxArray * superpixel_label_map_mx = prhs[0];
  const mxArray * boundary_map_mx = prhs[1];
  const mxArray * cascade_config_mx = prhs[2];
  
  const int * sizeImage;
  sizeImage = mxGetDimensions(superpixel_label_map_mx);
  int width = sizeImage[1], height = sizeImage[0];
  int n_pixel = height*width;
  unsigned int * superpixel_label_map_i = (unsigned int *)
    mxGetPr(superpixel_label_map_mx);
  plhs[0] = mxCreateNumericMatrix(height, width, mxUINT32_CLASS, mxREAL);
  unsigned int * superpixel_label_map = (unsigned int *) mxGetPr(plhs[0]);
  {
    int i;
    for(i=0; i<width*height; i++)
      superpixel_label_map[i] = superpixel_label_map_i[i];
  }
	
  BOUNDARY_VALUE * boundary_map = (BOUNDARY_VALUE *) mxGetPr(boundary_map_mx);
  
  // apply cascades
  Label_Pair_Id_Hash to_merge_pairs;

  Label_Pair_Id_Edge_Prop_Hash edge_prop;
  int n_stages = mxGetM(cascade_config_mx);
  mexPrintf("n_stages: %d\n", n_stages);
  int i;
  for(i=0; i<n_stages; i++){
    mxArray * stage_type_mx = mxGetCell(cascade_config_mx, i);
    mxArray * stage_param_mx = mxGetCell(cascade_config_mx, i+n_stages);
    char * stage_type = mxArrayToString(stage_type_mx);
    char * stage_param = mxArrayToString(stage_param_mx);

    mexPrintf("stage: %d\nstage_type: \"%s\"\nstage_param: \"%s\"\n",
              i, stage_type, stage_param);

    if(!strcmp(stage_type, "num_b<"))
      num_boundary_pixels_b_less_than(stage_type, stage_param,
                                      superpixel_label_map, boundary_map,
                                      width, height, edge_prop,
				      to_merge_pairs, i==0);
    if(!strcmp(stage_type, "num_b>"))
      num_boundary_pixels_b_greater_than(stage_type, stage_param,
                                      superpixel_label_map, boundary_map,
                                      width, height, edge_prop,
				      to_merge_pairs, i==0);
    
    mxFree(stage_type);
    mxFree(stage_param);
  }

  perform_mergers(superpixel_label_map, width*height, to_merge_pairs);

  mexPrintf("STOP: segment_cascade\n");
  return;
}

void perform_mergers(unsigned int * label_map, int n_pixel, Label_Pair_Id_Hash & to_merge_pairs){
  Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi> M(NULL);
  {
    int i;
    for(i=0; i<n_pixel; i++)
      M.add_new_set_inc(label_map[i]);
  }
  {
    Label_Pair_Id_Hash::iterator it;
    for(it=to_merge_pairs.begin(); it!=to_merge_pairs.end(); it++){
      if((*it).second==0)
	continue;
      Label_Pair p = id_2_label_pair((*it).first);
      M.merge(p.first, p.second);
    }
  }
  {
    int i;
    for(i=0; i<n_pixel; i++)
      label_map[i] = M.get_adam(label_map[i]);
  }
}

void num_boundary_pixels_b_less_than(char * stage_type, char * stage_param,
                                     unsigned int * superpixel_label_map,
                                     BOUNDARY_VALUE * boundary_map,
                                     int width, int height,
                                     Label_Pair_Id_Edge_Prop_Hash & edge_prop,
				     Label_Pair_Id_Hash & to_merge_pairs,
				     bool initialize_to_merge_pairs){

  edge_prop.clear();
  BOUNDARY_VALUE boundary_value_threshold;
  int n_pixel_threshold;
  int n_pixel_threshold_orient;
  {
    std::string sp(stage_param);
    std::istringstream spb(sp);
    int temp;
    spb >> temp;
    boundary_value_threshold = (BOUNDARY_VALUE) temp;
    spb >> n_pixel_threshold_orient;
    spb >> n_pixel_threshold;
    mexPrintf("boundary_value_threshold: %d\nn_pixel_threshold_orient: %d\nn_pixel_threshold: %d\n",
              (int) boundary_value_threshold, n_pixel_threshold_orient,
              n_pixel_threshold);
  }

  {
    int x, y;
    unsigned int * s;
    BOUNDARY_VALUE * b;
    Label_Pair_Id lp;
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
      }
    }
    for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
      }
    }
    for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
      }
    }
    for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
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
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;

          lp = label_pair_2_id(*(s-1), *(s+1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
	  /*
          lp = label_pair_2_id(*(s-height-1), *(s+height+1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;

          lp = label_pair_2_id(*(s-height+1), *(s+height-1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) < boundary_value_threshold;
	  */
        }
      }
      s++; // skip last element in column.
      b++;
    }
  }
  
  // merge (AND operation) with current to_merge_pairs
  {
    Label_Pair_Id_Edge_Prop_Hash::iterator ep_it;
    for(ep_it=edge_prop.begin(); ep_it!=edge_prop.end(); ep_it++){
      Label_Pair p = id_2_label_pair((*ep_it).first);
      if(p.first==p.second || p.first==0 || p.second==0)
        continue;
#ifdef __DEBUG__
      {
        mexPrintf("edge_prop[<%d,%d>].feature[0]: %d\n", p.first, p.second,
                  (*ep_it).second.boundary_feature[0]);
      }
#endif //__DEBUG__
      if(n_pixel_threshold_orient==0){
        if((*ep_it).second.boundary_feature[0]<n_pixel_threshold){
	  if(initialize_to_merge_pairs)
	    to_merge_pairs[(*ep_it).first]=1;
	}
	else
          to_merge_pairs[(*ep_it).first]=0;
      }
      else{
        if((*ep_it).second.boundary_feature[0]>n_pixel_threshold){
	  if(initialize_to_merge_pairs)
	    to_merge_pairs[(*ep_it).first]=1;
	}
	else
          to_merge_pairs[(*ep_it).first]=0;
      }
    }
#ifdef __DEBUG__
    {
      Label_Pair_Id_Hash::iterator it;
      mexPrintf("to_merge_pairs:\n");
      for(it=to_merge_pairs.begin(); it!=to_merge_pairs.end(); it++){
	if((*it).second==0)
	  continue;
	Label_Pair p = id_2_label_pair((*it).first);
        mexPrintf("%d %d\n", p.first, p.second);
      }
      mexPrintf("--\n");
    }
#endif //__DEBUG__
  }

}

void num_boundary_pixels_b_greater_than(char * stage_type, char * stage_param,
                                     unsigned int * superpixel_label_map,
                                     BOUNDARY_VALUE * boundary_map,
                                     int width, int height,
                                     Label_Pair_Id_Edge_Prop_Hash & edge_prop,
				     Label_Pair_Id_Hash & to_merge_pairs,
				     bool initialize_to_merge_pairs){
  mexPrintf("START: num_boundary_pixels_b_greater_than\n");
  edge_prop.clear();
  BOUNDARY_VALUE boundary_value_threshold;
  int n_pixel_threshold;
  int n_pixel_threshold_orient;
  {
    std::string sp(stage_param);
    std::istringstream spb(sp);
    int temp;
    spb >> temp;
    boundary_value_threshold = (BOUNDARY_VALUE) temp;
    spb >> n_pixel_threshold_orient;
    spb >> n_pixel_threshold;
    mexPrintf("boundary_value_threshold: %d\nn_pixel_threshold_orient: %d\nn_pixel_threshold: %d\n",
              (int) boundary_value_threshold, n_pixel_threshold_orient,
              n_pixel_threshold);
  }

  {
    int x, y;
    unsigned int * s;
    BOUNDARY_VALUE * b;
    Label_Pair_Id lp;
    for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;
      }
    }
    for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
        x<width-1;
        x++, s+=height, b+=height){
      if(*s==0){
        lp = label_pair_2_id(*(s-height), *(s+height));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;
      }
    }
    for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;
	if(lp==label_pair_2_id(2738, 6418) && (*b)>boundary_value_threshold)
	  mexPrintf("x:%d y:%d\n", 0, y);
      }
    }
    for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
        y<height-1;
        y++, s++, b++){
      if(*s==0){
        lp = label_pair_2_id(*(s-1), *(s+1));
        if(edge_prop.find(lp)==edge_prop.end()){
          edge_prop[lp].sum_boundary_value = 0;
          edge_prop[lp].n_boundary_pixel = 0;
          edge_prop[lp].boundary_feature[0] = 0;
        }
        edge_prop[lp].sum_boundary_value += *b;
        edge_prop[lp].n_boundary_pixel ++;
        edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;
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
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;

          lp = label_pair_2_id(*(s-1), *(s+1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold?1:0;
	  /*
          lp = label_pair_2_id(*(s-height-1), *(s+height+1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold;

          lp = label_pair_2_id(*(s-height+1), *(s+height-1));
          if(edge_prop.find(lp)==edge_prop.end()){
            edge_prop[lp].sum_boundary_value = 0;
            edge_prop[lp].n_boundary_pixel = 0;
            edge_prop[lp].boundary_feature[0] = 0;
          }
          edge_prop[lp].sum_boundary_value += *b;
          edge_prop[lp].n_boundary_pixel ++;
          edge_prop[lp].boundary_feature[0] += (*b) > boundary_value_threshold;
	  */
        }
      }
      s++; // skip last element in column.
      b++;
    }
  }
  
  // merge (AND operation) with current to_merge_pairs
  {
    Label_Pair_Id_Edge_Prop_Hash::iterator ep_it;
    for(ep_it=edge_prop.begin(); ep_it!=edge_prop.end(); ep_it++){
      Label_Pair p = id_2_label_pair((*ep_it).first);
      if(p.first==p.second || p.first==0 || p.second==0)
        continue;
#ifdef __DEBUG__
      {
        mexPrintf("edge_prop[<%d,%d>].feature[0]: %d\n", p.first, p.second,
                  (*ep_it).second.boundary_feature[0]);
      }
#endif //__DEBUG__
      if(n_pixel_threshold_orient==0){
        if((*ep_it).second.boundary_feature[0]<n_pixel_threshold){
	  if(initialize_to_merge_pairs)
	    to_merge_pairs[(*ep_it).first]=1;
	}
	else
          to_merge_pairs[(*ep_it).first]=0;
      }
      else{
        if((*ep_it).second.boundary_feature[0]>n_pixel_threshold){
	  if(initialize_to_merge_pairs)
	    to_merge_pairs[(*ep_it).first]=1;
	}
	else
          to_merge_pairs[(*ep_it).first]=0;
      }
    }
#ifdef __DEBUG__
    {
      Label_Pair_Id_Hash::iterator it;
      mexPrintf("to_merge_pairs:\n");
      for(it=to_merge_pairs.begin(); it!=to_merge_pairs.end(); it++){
	if((*it).second==0)
	  continue;
	Label_Pair p = id_2_label_pair((*it).first);
        mexPrintf("%d %d\n", p.first, p.second);
      }
      mexPrintf("--\n");
    }
#endif //__DEBUG__
  }
  mexPrintf("STOP: num_boundary_pixels_b_greater_than\n");

}
