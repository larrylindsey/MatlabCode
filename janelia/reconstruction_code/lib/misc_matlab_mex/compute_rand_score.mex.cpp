// compute rand score for to label mappings
// Shiv N. Vitaladevuni
//

#include <mex.h>
#include <math.h>
#include <ext/hash_map>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

namespace std{
  using namespace __gnu_cxx;
}

struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
		mexPrintf("\t1.	N0 x 2 uint32 matrix of mapping of N0 labels to M0 sets");
		mexPrintf("\t2.	N1 x 2 uint32 matrix of mapping of N1 labels to M1 sets");
    mexPrintf("output:\n");
		mexPrintf("\tRand score between the two mappings.\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

	const mxArray * map_0_mx = prhs[0];
	const mxArray * map_1_mx = prhs[1];
  
  int n_map_0 = mxGetM(map_0_mx);
  int n_map_1 = mxGetM(map_1_mx);

  Hash_UInt32_UInt32 hash_map_0, hash_map_1;
  {
    unsigned int i;
    unsigned int * map_0 = (unsigned int *) mxGetPr(map_0_mx);
    for(i=0; i<n_map_0; i++)
      hash_map_0[map_0[i]] = map_0[i+n_map_0];
    unsigned int * map_1 = (unsigned int *) mxGetPr(map_1_mx);
    for(i=0; i<n_map_1; i++)
      hash_map_1[map_1[i]] = map_1[i+n_map_1];
  }

  unsigned int rand_score_eq = 0;
  unsigned int rand_score_eq_bound = 0;
  unsigned int rand_score_neq = 0;
  unsigned int rand_score_neq_bound = 0;
  {
    Hash_UInt32_UInt32::iterator lp1, lp2;
    for(lp1=hash_map_0.begin(); lp1!=hash_map_0.end(); lp1++){
      lp2=lp1;
      lp2++;
      for(; lp2!=hash_map_0.end(); lp2++){
        bool eq_0 = (*lp1).second ==(*lp2).second;
        bool eq_1 = hash_map_1[(*lp1).first]==hash_map_1[(*lp2).first];
        rand_score_eq += eq_0 && !eq_1;
        rand_score_eq_bound += eq_0;
        rand_score_neq += !eq_0 && eq_1;
        rand_score_neq_bound += !eq_0;
      }
    }
  }

  plhs[0] = mxCreateDoubleMatrix(1, 4, mxREAL);
  *(mxGetPr(plhs[0])) = (double) rand_score_eq;
  *(mxGetPr(plhs[0])+1) = (double) rand_score_eq_bound;
  *(mxGetPr(plhs[0])+2) = (double) rand_score_neq;
  *(mxGetPr(plhs[0])+3) = (double) rand_score_neq_bound;
  return;
}
