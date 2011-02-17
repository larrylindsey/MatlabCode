// Get correspondence between segments  of two segmentation maps given the
// correspondences between their boundary elements
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//
// v0   08282008  init code
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*(p))
#define MAX_2_2(p,H) MAX(MAX(*(p),*(p+1)),MAX(*(p+H),*(p+H+1))) 

class Label_Pair{
	public:
		int label_0, label_1;
};
bool operator<(const Label_Pair& a, const Label_Pair& b){
	return (a.label_0<b.label_0)||((a.label_0==b.label_0)&&(a.label_1<b.label_1));
}

class Label_Pair_Stats{
	public:
		int label_0, label_1;
    int n_vote;
};

void get_two_labels(int index, int width, int height, double * label_map, int * label_a, int * label_b);

/*
 * input params:
 * 	1. watershed label map MxN matrix
 *	2. boundary map
 *	3. threshold on mean boundary value below which segments are merged.
 *
 * output:
 *	1. merged watershed label map
 *	2. mapping from input segment labels to merged segment labels
 */
void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage seg0_to_seg1_map =\n");
    mexPrintf("\tcorrespond_segments_from_boundary_correspondence(boundary0_to_boundary1_map, ...\n");
    mexPrintf("\t\tsegment_map_0, segment_map_1)\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. boundary0_to_boundary1_map: Map of boundary indexes of\n");
    mexPrintf("\t\tsegment_map_0 to segment_map_1 [2xN] matrix. Indexes as generated\n");
    mexPrintf("\t\tby sub2ind.\n");
		mexPrintf("\t2. segmentation_map_0 label map with 0's for boundaries - PxQ matrix.\n");
		mexPrintf("\t3. segmentation_map_1 PxQ matrix\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. mapping from labels from segment0 to segment1. Has 0's when no correspondence.\n");
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

	int numDim = mxGetNumberOfDimensions(prhs[1]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 2\n");
		return;
	}

	double * b_2_b_map = mxGetPr(prhs[0]);
	double * segment_map_0 = mxGetPr(prhs[1]);
	double * segment_map_1 = mxGetPr(prhs[2]);
  
	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[1]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	
  // find the number of pairs of pixels from segment0 and segment1 with nonzero labels
  // give the maximum number of pairs of labels that can correspond
	int max_n_segment_overlap_pair = 0;
  double * s_m_0 = segment_map_0, * s_m_1 = segment_map_1;
	for(int i=0; i<n_pixel; i++, s_m_0++, s_m_1++)
		if(*s_m_0>0 && *s_m_1>0)
      max_n_segment_overlap_pair++;
  
  // collect segment labels pairs for matched boundary elements
  Label_Pair * label_pairs = (Label_Pair *) new Label_Pair[max_n_segment_overlap_pair];
  Label_Pair * l_p = label_pairs;

  /* for debugging
  double a[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int b, c;
  get_two_labels(4, 3, 3, a, &b, &c);
  mexPrintf("%d %d\n", b, c);
  return;
  */
  
  const int * size_temp = mxGetDimensions(prhs[0]);
  int n_matched_boundary = size_temp[1];
  int index_0, index_1;
  int label_0a, label_0b; // two segments for a boundary
  int label_0t; // segment label in the other label map
  int label_1a, label_1b; // two segments for a boundary
  int label_1t; // segment label in the other label map
  int n_segment_overlap_pair = 0;
  for(int i=0; i<2*n_matched_boundary; i+=2){
    index_0 = (int)b_2_b_map[i]-1; // 0 shifted
    index_1 = (int)b_2_b_map[i+1]-1;

    label_0t = (int)*(segment_map_1+index_0);
    label_1t = (int)*(segment_map_0+index_1);
    if(label_0t>0 && label_1t>0){
      get_two_labels(index_0, width, height, segment_map_0, & label_0a, & label_0b);
      get_two_labels(index_1, width, height, segment_map_1, & label_1a, & label_1b);
      if(label_0a<=0 || label_0b<=0 || label_1a<=0 || label_1b<=0)
        continue;
      if(label_0t!=label_1a && label_0t!=label_1b)
        continue;
      if(label_1t!=label_0a && label_1t!=label_0b)
        continue;
      
      if(label_1t==label_0a && label_0t==label_1a){
        l_p->label_0 = label_0b;
        l_p->label_1 = label_1a;
        l_p++;
        n_segment_overlap_pair++;
        l_p->label_0 = label_0a;
        l_p->label_1 = label_1b;
        l_p++;
        n_segment_overlap_pair++;
      }
      if(label_1t==label_0b && label_0t==label_1a){
        l_p->label_0 = label_0a;
        l_p->label_1 = label_1a;
        l_p++;
        n_segment_overlap_pair++;
        l_p->label_0 = label_0b;
        l_p->label_1 = label_1b;
        l_p++;
        n_segment_overlap_pair++;
      }
      if(label_1t==label_0a && label_0t==label_1b){
        l_p->label_0 = label_0b;
        l_p->label_1 = label_1b;
        l_p++;
        n_segment_overlap_pair++;
        l_p->label_0 = label_0a;
        l_p->label_1 = label_1a;
        l_p++;
        n_segment_overlap_pair++;
      }
      if(label_1t==label_0b && label_0t==label_1b){
        l_p->label_0 = label_0a;
        l_p->label_1 = label_1b;
        l_p++;
        n_segment_overlap_pair++;
        l_p->label_0 = label_0b;
        l_p->label_1 = label_1a;
        l_p++;
        n_segment_overlap_pair++;
      }
    }
    else{
      if(index_0==index_1){
      }
    }
  }
  
  Label_Pair * sorted_label_pairs = label_pairs;
	std::sort(sorted_label_pairs, sorted_label_pairs+n_segment_overlap_pair-1);
  
  // get statistics for each label pair
  Label_Pair_Stats * label_pair_stats = (Label_Pair_Stats *) new Label_Pair_Stats[n_segment_overlap_pair];
  Label_Pair_Stats * l_p_s = label_pair_stats;
  l_p = sorted_label_pairs;
	int prev_label_0=0, prev_label_1=0;
  int n_vote=0;
  int n_label_pair=0;
  for(int i=0; i<n_segment_overlap_pair; i++, l_p++){
		if(l_p->label_0<=0 || l_p->label_1<=0)
			continue;
		
		if(l_p->label_1!=prev_label_1 || l_p->label_0!=prev_label_0){
      if(n_vote!=0){
        l_p_s->n_vote = n_vote;
        l_p_s++;
      }
      l_p_s->label_0 = l_p->label_0;
      l_p_s->label_1 = l_p->label_1;
      n_vote=1;
      n_label_pair++;
			
			prev_label_0=l_p->label_0;
			prev_label_1=l_p->label_1;
		}
    else{
      n_vote++;
    }
	}
  if(n_vote!=0){
    l_p_s->n_vote = n_vote;
  }
  
  #define N_OUT_DIM 3
	plhs[0] = mxCreateDoubleMatrix(N_OUT_DIM, n_label_pair, mxREAL);
  double * pairwise_votes = mxGetPr(plhs[0]);
  l_p_s = label_pair_stats;
  for(int i=0; i<n_label_pair*N_OUT_DIM; i+=N_OUT_DIM, l_p_s++){
    pairwise_votes[i] = l_p_s->label_0;
    pairwise_votes[i+1] = l_p_s->label_1;
    pairwise_votes[i+2] = l_p_s->n_vote;
  }
	
  delete [] label_pairs;
  
	return;
}

void get_two_labels(int index, int width, int height, double * label_map, int * label_a, int * label_b){
  *label_a=-10;
  *label_b=-10;
  int x = index / height;
  int y = index % height;
  int * l = label_a;
  label_map += index;
  if(x>0){
    if(y>0)
      if((*(label_map-1-height))>0){
        *l = (int)*(label_map-1-height);
        l = label_b;
      }
    if((*(label_map-height))>0 && *label_a!=*(label_map-height)){
      *l = (int)*(label_map-height);
      l = label_b;
    }
    if(y<height-1)
      if((*(label_map+1-height))>0 && *label_a!=*(label_map+1-height)){
        *l = (int)*(label_map+1-height);
        l = label_b;
      }
  }
  
  if(y>0)
    if((*(label_map-1))>0 && *label_a!=*(label_map-1)){
      *l = (int)*(label_map-1);
      l = label_b;
    }
  if((*label_map)>0 && *label_a!=*label_map){
    *l = (int)*label_map;
    l = label_b;
  }
  if(y<height-1)
    if((*(label_map+1))>0 && *label_a!=*(label_map+1)){
      *l = (int)*(label_map+1);
      l = label_b;
    }
  
  if(x<width-1){
    if(y>0)
      if((*(label_map-1+height))>0 && *label_a!=*(label_map-1+height)){
        *l = (int)*(label_map-1+height);
        l = label_b;
      }
    if((*(label_map+height))>0 && *label_a!=*(label_map+height)){
      *l = (int)*(label_map+height);
      l = label_b;
    }
    if(y<height-1)
      if((*(label_map+1+height))>0 && *label_a!=*(label_map+1+height)){
        *l = (int)*(label_map+1+height);
        l = label_b;
      }
  }
  return;
}
