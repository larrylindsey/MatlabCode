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

#include <stdio.h>
#include "../lib/maxflow-v3.0.src/graph.cpp"
#include "../lib/maxflow-v3.0.src/maxflow.cpp"

void err_f(char * msg){
  mexPrintf("%s", msg);
}

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
		
  typedef Graph<int,int,int> GraphType;
	GraphType *g = new GraphType(n_pixel, n_pixel*2, err_f); 

  for(int i=0; i<n_pixel; i++)
  	g -> add_node(); 

  int * prior_cost = (int *) mxGetPr(label_energy_prior_mx);
  for(int i=0; i<n_pixel; i++)
    g -> add_tweights( i,  prior_cost[2*i], prior_cost[2*i+1]);
  
  int * b_m_p = (int *) mxGetPr(boundary_map_mx);
  int node_id=0;
  int weight;
  for(int x=0; x<width; x++){
    for(int y=0; y<height-1; y++){
      weight = ((*b_m_p)+(*(b_m_p+1)))/2;
      g -> add_edge(node_id, node_id+1, weight, weight);
      node_id++;
      b_m_p++;
    }
    node_id++;
    b_m_p++;
  }
      
  b_m_p = (int *) mxGetPr(boundary_map_mx);
  node_id=0;
  for(int x=0; x<width-1; x++){
    for(int y=0; y<height; y++){
      weight = ((*b_m_p)+(*(b_m_p+height)))/2;
      g -> add_edge(node_id, node_id+height, weight, weight);
      node_id++;
      b_m_p++;
    }
  }

	int flow = g -> maxflow();

  plhs[0] = mxCreateDoubleMatrix(height, width, mxREAL);
  double * l = mxGetPr(plhs[0]);
  for(int i=0; i<n_pixel; i++, l++)
    *l = g->what_segment(i)==GraphType::SOURCE;

	delete g;

	return;
}
