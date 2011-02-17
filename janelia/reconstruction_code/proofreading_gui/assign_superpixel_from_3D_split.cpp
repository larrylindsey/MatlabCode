// Reassign superpixels based on a 3D split
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	11172008	init. code
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
		unsigned short int superpixel_id;
    int body_id;
};
bool operator<(const Label_Pair& a, const Label_Pair& b){
	return (a.superpixel_id<b.superpixel_id)||((a.superpixel_id==b.superpixel_id)&&(a.body_id<b.body_id));
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage superpixel_2_body_map = assign_superpixel_from_3D_split(cat, label_new, ...\n\timage_stack_scale, z);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. cat{}: cell array of P 2D label maps volume sMxsN uint16 matrices\n");
		mexPrintf("\t2. label_new_stack: 3D stack of new label volume MxNxP int32 matrix\n");
		mexPrintf("\t3. image_stack_scale (s): down-scaling along x and y axis\n");
		mexPrintf("\t4. z: z plane to consider\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. mapping from superpixels to bodies 2xR\n");
		return;
	}
	if(nrhs!=4){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>2){
		mexErrMsgTxt("Wrong number of outputs\n");
		return;
	}

  const mxArray * cat_mx = prhs[0];
  const mxArray * label_new_stack_mx = prhs[1];
  const mxArray * image_stack_scale_mx = prhs[2];
  const mxArray * z_mx = prhs[3];

	const int * size_0 = mxGetDimensions(cat_mx);
  int n_plane = MAX(size_0[0], size_0[1]);
  const int * size_0_0 = mxGetDimensions(mxGetCell(cat_mx, 0));
  int original_width = size_0_0[0];
  
  int z = (int) * mxGetPr(z_mx);
  if(z<0 || z>=n_plane){
		mexErrMsgTxt("Invalid z\n");
		return;
  }

  unsigned short int * superpixel_map = (unsigned short int *) mxGetPr(mxGetCell(cat_mx, z));
  
  int * label_new = (int *) mxGetPr(label_new_stack_mx);
	const int * size_1 = mxGetDimensions(label_new_stack_mx);
	const int width = size_1[0], height = size_1[1], depth = size_1[2];
  
  int image_stack_scale = (int) * mxGetPr(image_stack_scale_mx);
  
  // Get votes of superpixel to body assignments
  unsigned short int sp_id;
  int * l_n = label_new + width*height*z;
  int n_pair = 0;
  Label_Pair * label_pairs = (Label_Pair *) new Label_Pair[width*height];
  Label_Pair * l_p = label_pairs;
  for(int i=0; i<width*height; i++, l_n++){
    if(*l_n<=0)
      continue;
    sp_id = superpixel_map[image_stack_scale*(original_width*(i/width) + (i%width))];
    if(sp_id==0)
      continue;
    
    l_p->superpixel_id = sp_id;
    l_p->body_id = *l_n;
    l_p++;
    n_pair++;
  }
  
  if(n_pair==0){
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    return;
  }
  
	// Sort the vote pairs
	std::sort(label_pairs, label_pairs+n_pair);
	
	l_p = label_pairs;
	int prev_superpixel_id=0;
  int n_superpixel = 0;
	for(int i=0; i<n_pair; i++, l_p++)
		if(l_p->superpixel_id!=prev_superpixel_id)
      n_superpixel++;
  
  plhs[0] = mxCreateDoubleMatrix(2, n_superpixel, mxREAL);
  double * superpixel_to_body_map = mxGetPr(plhs[0]);
  
  double * sp_2_b_m = superpixel_to_body_map;
  l_p = label_pairs;
  prev_superpixel_id = 0;
  Label_Pair * l_p_first_element = NULL;
  int n_vote = 0;
	for(int i=0; i<n_pair; i++, l_p++){
		if(l_p->superpixel_id!=prev_superpixel_id){
      if(n_vote!=0){
        *(sp_2_b_m) = (l_p_first_element+(int)(n_vote/2))->superpixel_id;
        *(sp_2_b_m+1) = (l_p_first_element+(int)(n_vote/2))->body_id;
        sp_2_b_m += 2;
      }
      
      n_vote=1;
      l_p_first_element = l_p;
			prev_superpixel_id = l_p->superpixel_id;
		}
    else{
      n_vote++;
    }
	}
  *(sp_2_b_m) = (l_p_first_element+(int)(n_vote/2))->superpixel_id;
  *(sp_2_b_m+1) = (l_p_first_element+(int)(n_vote/2))->body_id;
	
	delete [] label_pairs;
		
	return;
}


