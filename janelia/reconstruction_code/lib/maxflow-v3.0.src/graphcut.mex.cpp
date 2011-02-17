// Wrapper for Graph Cut
//
// For Graph Cut citation and licensing see $(LIB_DIR)/MRF2.1/README.txt
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
#include "graph.cpp"
#include "maxflow.cpp"

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
		mexPrintf("Usage labels = graphcut(label_prior, pairwise_smooth, n_iteration)\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. label priors - Mx2 <cost for label 0, cost for label 1> (int32)\n");
		mexPrintf("\t2. smoothing term - Rx3 <node id 1, node id 2, smoothing weight> (int32)'\n");
		mexPrintf("\t3. number of iterations\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. labels\n");
		return;
	}
	if(nrhs!=3){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}

  const mxArray * label_prior_mx = prhs[0];
  const mxArray * pairwise_smooth_mx = prhs[1];
  const int n_iteration = (int) *mxGetPr(prhs[2]);
  
	if(!mxIsInt32(label_prior_mx)){
		mexErrMsgTxt("Label prior must by Int32 values\n");
		return;
	}
	if(!mxIsInt32(pairwise_smooth_mx)){
		mexErrMsgTxt("Pairwise smoothing terms must by Int32 values\n");
		return;
	}

  const int n_node = mxGetM(label_prior_mx);
  mexPrintf("n_node:%d\n", n_node);
	const int n_pair = mxGetM(pairwise_smooth_mx);
		
  typedef Graph<int,int,int> GraphType;
	GraphType *g = new GraphType(n_node, n_pair, err_f); 

  for(int i=0; i<n_node; i++)
  	g -> add_node(); 

  int * label_prior = (int *) mxGetPr(label_prior_mx);
  for(int i=0; i<n_node; i++)
    g -> add_tweights( i,  label_prior[i], label_prior[i+n_node]);
  
  int * pairwise_smooth = (int *) mxGetPr(pairwise_smooth_mx);
  for(int p=0; p<n_pair; p++){
    g -> add_edge(pairwise_smooth[p], pairwise_smooth[p+n_pair],
                  pairwise_smooth[p+n_pair*2], pairwise_smooth[p+n_pair*2]);
  }
      
	int flow = g -> maxflow();

  plhs[0] = mxCreateDoubleMatrix(n_node, 1, mxREAL);
  double * l = mxGetPr(plhs[0]);
  for(int i=0; i<n_node; i++, l++)
    *l = g->what_segment(i)==GraphType::SOURCE;

	delete g;

	return;
}
