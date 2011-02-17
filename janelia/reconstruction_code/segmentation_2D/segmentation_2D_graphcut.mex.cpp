// compute segmentation on a boundary map using Graph Cut given seeds of for label 0 and 1
//
// For Graph Cut citation see $(LIB_DIR)/MRF2.1/README.txt
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	04012008	init. code
//

#include <mex.h>
#include <math.h>
#include <algorithm>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

#define MAX_3_3(p,H) MAX(MAX(MAX(MAX(*(p-1-H),*(p-1)),MAX(*(p-1+H),*(p+H))),MAX(MAX(*(p+1+H),*(p+1)),MAX(*(p+1-H),*(p-H)))),*p)
#define MAX_2_2(p,H) MAX(MAX(*p,*(p+1)),MAX(*(p+H),*(p+H+1))) 

#include <mrf.h>
#include <GCoptimization.h>

MRF::CostVal * prior_cost;
MRF::CostVal * smooth_cost;
MRF::CostVal * edge_weight_horz;
MRF::CostVal * edge_weight_vert;
DataCost * data;
EnergyFunction * construct_energy_function(const mxArray * boundary_map_mx, const mxArray * label_energy_prior_mx);

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage label_map = segmentation_2D_graphcut(boundary_map, label_energy_prior)\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. boundary field on which graphcut is to be performed - MxN matrix\n");
		mexPrintf("\t2. prior label energy for each pixel - integer array 2xMN matrix <prob. label 0, prob. label 1>'\n");
		mexPrintf("\t3. number of iterations\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. label map\n");
		return;
	}
	if(nrhs!=3){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}

  const mxArray * boundary_map_mx = prhs[0];
  const mxArray * label_energy_prior_mx = prhs[1];
  const int n_iteration = (int) *mxGetPr(prhs[2]);
  
  if(sizeof(MRF::CostVal)!=sizeof(int)){
    mexErrMsgTxt("Graph Cut's CostVal type is different from Int32 .. exiting\n");
		return;
	}
  
	int numDim = mxGetNumberOfDimensions(boundary_map_mx);
	if(numDim>2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
		return;
	}

	if(!mxIsInt32(boundary_map_mx)){
		mexErrMsgTxt("Boundary costs must by Int32 values\n");
		return;
	}
	if(!mxIsInt32(label_energy_prior_mx)){
		mexErrMsgTxt("Label energies must by Int32 values\n");
		return;
	}
  
  const int * sizeImage;
	sizeImage = mxGetDimensions(boundary_map_mx);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
		
  MRF* mrf;
  EnergyFunction *energy = construct_energy_function(boundary_map_mx, label_energy_prior_mx);
  MRF::EnergyVal E;

  mrf = new Swap(height,width,2,energy);
  mrf->initialize();
  mrf->clearAnswer();
  
  E = mrf->totalEnergy();
  //printf("Energy at the Start= %g (%g,%g)\n", (float)E, (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());
  
  #ifdef COUNT_TRUNCATIONS
  truncCnt = totalCnt = 0;
  #endif
  float tot_t = 0;
  float t;
  for(int iter=0; iter<n_iteration; iter++) {
    mrf->optimize(1, t);
    
    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    //printf("energy = %g (%f secs)\n", (float)E, tot_t);
  }
  
  MRF::Label * label = mrf->getAnswerPtr();
  plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
  double * l_m_p = mxGetPr(plhs[0]);
  for(int i=0; i<n_pixel; i++, l_m_p++,label++)
    *l_m_p = (double)(*label);
  
  delete mrf;
  delete energy;
  delete data;
  delete [] edge_weight_horz;
  delete [] edge_weight_vert;
  delete [] smooth_cost;
	return;
}


EnergyFunction * construct_energy_function(const mxArray * boundary_map_mx, const mxArray * label_energy_prior_mx){
  
	const int * sizeImage;
	sizeImage = mxGetDimensions(boundary_map_mx);
	int width = sizeImage[1], height = sizeImage[0];
  
	int n_pixel = height*width;
  
  // single-variable terms 
  prior_cost = (MRF::CostVal *) mxGetPr(label_energy_prior_mx);
  data = new DataCost(prior_cost);
  
  // pairwise terms
  smooth_cost = (MRF::CostVal *) new MRF::CostVal[4];
  smooth_cost[0] = 0; smooth_cost[1] = 1; smooth_cost[2] = 1; smooth_cost[3] = 0;
  
  edge_weight_horz = (MRF::CostVal *) new MRF::CostVal[n_pixel];
  edge_weight_vert = (MRF::CostVal *) new MRF::CostVal[n_pixel];

  MRF::CostVal * e_w_p = edge_weight_horz;
  int * b_m_p = (int *) mxGetPr(boundary_map_mx);
  for(int x=0; x<width; x++, e_w_p++, b_m_p++){
    for(int y=0; y<height-1; y++, e_w_p++, b_m_p++){
      *e_w_p = ((*b_m_p)+(*(b_m_p+1)))/2;
    }
  }
      
  e_w_p = edge_weight_vert;
  b_m_p = (int *) mxGetPr(boundary_map_mx);
  for(int x=0; x<width-1; x++){
    for(int y=0; y<height; y++, e_w_p++, b_m_p++){
      *e_w_p = ((*b_m_p)+(*(b_m_p+height)))/2;
    }
  }
  
  SmoothnessCost *smooth = (SmoothnessCost *) new SmoothnessCost(smooth_cost,edge_weight_horz,edge_weight_vert);
  
  EnergyFunction * energy = (EnergyFunction *) new EnergyFunction(data,smooth);
  
  return energy;
}

