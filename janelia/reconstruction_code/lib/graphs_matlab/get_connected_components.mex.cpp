// Get connected components in a weighted undirected graph given it's
// sparse adjacency matrix. A threshold is applied on the edge
// weights for them to be considered for adjacency.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//

#include <mex.h>
#include <queue>


#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
		mexPrintf("\t1.	MxM double sparse matrix\n");
		mexPrintf("\t2.	Threshold on the edge weight\n");
    mexPrintf("output:\n");
		mexPrintf("\t1. Mx1 uint32 Connected component ids\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  const mxArray * A_mx = prhs[0];
  const mxArray * threshold_mx = prhs[1];

  const double * A = mxGetPr(A_mx);
//  mwIndex * ir = mxGe+tIr(A_mx);
//  mwIndex * jc = mxGetJc(A_mx);
  int * ir = mxGetIr(A_mx);
  int * jc = mxGetJc(A_mx);

  const unsigned int n_node = (unsigned int) mxGetN(A_mx);
  mexPrintf("number of nodes: %d\n", n_node);
  const double threshold = * mxGetPr(threshold_mx);
  
  const int size_temp[2] = {n_node, 1};
  plhs[0] = mxCreateNumericArray(2, size_temp, mxUINT32_CLASS,
                                 mxREAL);
  unsigned int * component_id = (unsigned int *) mxGetPr(plhs[0]);

  std::queue<unsigned int> Q;
  unsigned int n, curr_component_id=1;
  unsigned int n_neigh, edge_id, neigh_id, i, node_id;
  for(n = 0; n<n_node; n++){
    if(component_id[n]>0)
      continue;
    component_id[n] = curr_component_id;
    for(edge_id=jc[n]; edge_id<jc[n+1]; edge_id++){
      if(A[edge_id]<=threshold)
        continue; // not considered for adjacency
      neigh_id = ir[edge_id];
      if(component_id[neigh_id]>0)
        continue; // already processed.
      Q.push(neigh_id);
    }
    while(!Q.empty()){
      node_id = Q.front();
      Q.pop();
      component_id[node_id] = curr_component_id;
      for(edge_id=jc[node_id]; edge_id<jc[node_id+1]; edge_id++){
        if(A[edge_id]<=threshold)
          continue; // not considered for adjacency
        neigh_id = ir[edge_id];
        if(component_id[neigh_id]>0)
          continue; // already processed.
        Q.push(neigh_id);
      }
    }
    curr_component_id++;
  }
  return;
}
