
#include <mex.h>
#include <algorithm>                     // for std::for_each
#include <ext/hash_map>
#include <boost/utility.hpp>             // for boost::tie

// CSA++
#include "Array.hh"
#include "kofn.hh"
#include "csa.hh"

namespace std{
  using namespace __gnu_cxx;
}

// Hash function from labels to unsigned int
typedef unsigned int Label;

struct equui{ bool operator()(const Label s1, const Label s2) const
    {return s1==s2;}};

typedef std::hash_map<Label, unsigned int, std::hash<Label>, equui> Label_Hash;

// Hash function from pair of labels to unsigned int
typedef unsigned long int Label_Pair;

Label_Pair label_pair_2_id(Label s1, Label s2){
  Label_Pair lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair lp;
  lp = s2;
  lp |= lp1;
  return lp;
}

std::pair<Label, Label> id_2_label_pair(Label_Pair lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

struct equli{ bool operator()(const Label_Pair s1, const Label_Pair s2) const
    {return s1==s2;}};

typedef std::hash_map<Label_Pair, unsigned int, std::hash<Label_Pair>, equli> Label_Pair_Hash;


// CSA code needs integer weights.  Use this multiplier to convert
// floating-point weights to integers.
static const int multiplier = 100;

// The degree of outlier connections.
static const int degree = 300;

// Outlier cost
static int outlierCost = 100000;

static const int  WEIGHT_SCALE = 10;

using namespace boost;

void assert1(bool query, int error_id){
  if(!query){
    char msg[256];
    sprintf(msg, "Error [%d]\n", error_id);
   	mexErrMsgTxt(msg);
  }
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage diff_stack = compare_label_stacks(stack_0, stack_1);\n");
		return;
	}
	if(nrhs!=3){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs!=1 && nlhs!=2){
		mexErrMsgTxt("Wrong number of outputs\n");
		return;
	}

	const int * size_0 = mxGetDimensions(prhs[0]);
	int depth = size_0[2], height = size_0[0], width = size_0[1];
	const int * size_1 = mxGetDimensions(prhs[1]);
	int depth_1 = size_1[2], height_1 = size_1[0], width_1 = size_1[1];
  if(depth!=depth_1 || height!=height_1 || width!=width_1){
    mexErrMsgTxt("Stack dimensions don't match\n");
    return;
  }

  double min_overlap_norm_threshold = * mxGetPr(prhs[2]);
  
  int n_pixel_plane = height * width;
  int n_pixel = n_pixel_plane * depth;
  mexPrintf("n_pixel_plane: %d, n_pixel: %d\n", n_pixel_plane, n_pixel);
  
	Label * stack_0 = (Label *) mxGetPr(prhs[0]);
	Label * stack_1 = (Label *) mxGetPr(prhs[1]);

  mexPrintf("get the number of nodes in A and B\n");
  Label_Hash label_2_node_id_0, label_2_node_id_1;
  std::vector<Label> node_id_2_label_0, node_id_2_label_1;
  Label_Hash n_pixel_label_0, n_pixel_label_1;
  {
    int i;
    for(i=0; i<n_pixel; i++){
      if(stack_0[i]==0)
        continue;
      if(n_pixel_label_0[stack_0[i]]==0){
        label_2_node_id_0[stack_0[i]] = node_id_2_label_0.size();
        node_id_2_label_0.push_back(stack_0[i]);
      }        
      n_pixel_label_0[stack_0[i]]++;
    }
    for(i=0; i<n_pixel; i++){
      if(stack_1[i]==0)
        continue;
      if(n_pixel_label_1[stack_1[i]]==0){
        label_2_node_id_1[stack_1[i]] = node_id_2_label_1.size();
        node_id_2_label_1.push_back(stack_1[i]);
      }
      n_pixel_label_1[stack_1[i]]++;
    }
  }
  int n_label_0, n_label_1;
  n_label_0 = n_pixel_label_0.size();
  n_label_1 = n_pixel_label_1.size();
  mexPrintf("n_label_0:%d n_label_1:%d\n", n_label_0, n_label_1);
  mexPrintf("node_id_2_label_0.size:%d node_id_2_label_1.size:%d\n",
            node_id_2_label_0.size(), node_id_2_label_1.size());
  
  unsigned int max_overlap = 0;
  Label_Pair_Hash n_pixel_label_pair;
  {
    mexPrintf("count label pairs - determine the weights in bipartite graph\n");
    int i;
    for(i=0; i<n_pixel; i++){
      if(stack_0[i]==0 || stack_1[i]==0)
        continue;
      n_pixel_label_pair[label_pair_2_id(stack_0[i], stack_1[i])]++;
    }

#ifdef __DEBUG__
    {
      mexPrintf("Printing label pairs and overlap areas:\n");
      Label_Pair_Hash::iterator lp_it;
      Label_Pair lp;
      Label l0, l1;
      for(lp_it=n_pixel_label_pair.begin(); lp_it!=n_pixel_label_pair.end();
          lp_it++){
        lp = (*lp_it).first;
        l1 = lp & 0X00000000ffffffff;
        l0 = lp >> 32;
        mexPrintf("%u %u %u\n", l0, l1, (*lp_it).second);
      }
      mexPrintf("done.\n");
    }
#endif

    mexPrintf("normalizing overlap areas for graph weights.\n");
    Label_Pair_Hash::iterator lp_it;
    for(lp_it=n_pixel_label_pair.begin(); lp_it!=n_pixel_label_pair.end();
        lp_it++)
      max_overlap = std::max(max_overlap, (*lp_it).second); 

    for(lp_it=n_pixel_label_pair.begin(); lp_it!=n_pixel_label_pair.end();
        lp_it++){
//      mexPrintf("%d ", (*lp_it).second);
      (*lp_it).second = std::max(1, (int)max_overlap - (int)(*lp_it).second);
//      (*lp_it).second /= WEIGHT_SCALE;
//      mexPrintf("%d\n", (*lp_it).second);
    }

    outlierCost = max_overlap;
    
#ifdef __DEBUG__
    {
      mexPrintf("Printing label pairs and bipartite graph weights:\n");
      Label_Pair_Hash::iterator lp_it;
      Label_Pair lp;
      Label l0, l1;
      for(lp_it=n_pixel_label_pair.begin(); lp_it!=n_pixel_label_pair.end();
          lp_it++){
        lp = (*lp_it).first;
        l1 = lp & 0X00000000ffffffff;
        l0 = lp >> 32;
        mexPrintf("%u %u %u\n", l0, l1, (*lp_it).second);
      }
      mexPrintf("done.\n");
    }
#endif
  }
  mexPrintf("n_pixel_label_pair.size:%d\n", n_pixel_label_pair.size());
  
  mexPrintf("Computing optimal assignment between segments.\n");
  std::vector<std::pair<std::pair<Label, Label>, unsigned int> > matched_pair;
  {
    // 
    // Code in this scope has been predominantly copied and modified
    // from segbench/Benchmark/match.cc. "segbench" code provided by
    // Jitendra Malik's group at UC Berkeley. Cite Martin et al. PAMI
    // 2004 if using any part of this code.
    //

    mexPrintf("Building bipartite graph\n");
    // The cardinality of the match is n.
    const int n = n_label_0 + n_label_1;
    const int nmin = std::min(n_label_0,n_label_1);
    const int nmax = std::max(n_label_0,n_label_1);

    // Compute the degree of various outlier connections.
    const int d1 = std::max(0,std::min(degree,n_label_0-1)); // from map1
    const int d2 = std::max(0,std::min(degree,n_label_1-1)); // from map2
    const int d3 = std::min(degree,std::min(n_label_0,n_label_1)); // between outliers
    const int dmax = std::max(d1,std::max(d2,d3));

    assert1 (n_label_0 == 0 || (d1 >= 0 && d1 < n_label_0), 1);
    assert1 (n_label_1 == 0 || (d2 >= 0 && d2 < n_label_1), 2);
    assert1 (d3 >= 0 && d3 <= nmin, 3);

    // Count the number of edges.
    int m = 0;
    m += n_pixel_label_pair.size(); 	// real connections
    m += d1 * n_label_0;	// outlier connections
    m += d2 * n_label_1;	// outlier connections
    m += d3 * nmax;	// outlier-outlier connections
    m += n; 		// high-cost perfect match overlay

    // If the graph is empty, then there's nothing to do.
    if (m == 0) {
      return;
    }

    // Weight of outlier connections.
    const int ow = outlierCost;// * multiplier;

    // Scratch array for outlier edges.
    Array1D<int> outliers (dmax);

    // Construct the input graph for the assignment problem.
    Array2D<int> igraph (m,3);
    int count = 0;
    // real edges
    Label_Pair_Hash::iterator lp_it;
    Label_Pair lp;
    std::pair<Label,Label> lps;
    for (lp_it=n_pixel_label_pair.begin(); lp_it!=n_pixel_label_pair.end();
         lp_it++) {
      lp = (*lp_it).first;
      lps = id_2_label_pair(lp);
      int i = label_2_node_id_0[lps.first];
      int j = label_2_node_id_1[lps.second];
//      mexPrintf("%d %d %d\n", i, j, (*lp_it).second);
      assert1 (i >= 0 && i < n_label_0, 4);
      assert1 (j >= 0 && j < n_label_1, 5);
      igraph(count,0) = i;
      igraph(count,1) = j;
      igraph(count,2) = (*lp_it).second;
      count++;
    }
    // outliers edges for map1, exclude diagonal
    for (int i = 0; i < n_label_0; i++) {
      kOfN(d1,n_label_0-1,outliers.data());
      for (int a = 0; a < d1; a++) {
        int j = outliers(a);
        if (j >= i) { j++; }
        assert1 (i != j, 6);
        assert1 (j >= 0 && j < n_label_0, 7);
        igraph(count,0) = i;
        igraph(count,1) = n_label_1 + j;
        igraph(count,2) = ow;
        count++;
      }
    }
    // outliers edges for map2, exclude diagonal
    for (int j = 0; j < n_label_1; j++) {
      kOfN(d2,n_label_1-1,outliers.data());
      for (int a = 0; a < d2; a++) {
        int i = outliers(a);
        if (i >= j) { i++; }
        assert1 (i != j, 8);
        assert1 (i >= 0 && i < n_label_1, 9);
        igraph(count,0) = n_label_0 + i;
        igraph(count,1) = j;
        igraph(count,2) = ow;
        count++;
      }
    }
    // outlier-to-outlier edges
    for (int i = 0; i < nmax; i++) {
      kOfN(d3,nmin,outliers.data());
      for (int a = 0; a < d3; a++) {
        const int j = outliers(a);
        assert1 (j >= 0 && j < nmin, 10);
        if (n_label_0 < n_label_1) {
          assert1 (i >= 0 && i < n_label_1, 11);
          assert1 (j >= 0 && j < n_label_0, 12);
          igraph(count,0) = n_label_0 + i;
          igraph(count,1) = n_label_1 + j;
        } else {
          assert1 (i >= 0 && i < n_label_0, 13);
          assert1 (j >= 0 && j < n_label_1, 14);
          igraph(count,0) = n_label_0 + j;
          igraph(count,1) = n_label_1 + i;
        }
        igraph(count,2) = ow;
        count++;
      }
    }
    // perfect match overlay (diagonal)
    for (int i = 0; i < n_label_0; i++) {
      igraph(count,0) = i;
      igraph(count,1) = n_label_1 + i;
      igraph(count,2) = ow * multiplier;
      count++;
    }
    for (int i = 0; i < n_label_1; i++) {
      igraph(count,0) = n_label_0 + i;
      igraph(count,1) = i;
      igraph(count,2) = ow * multiplier;
      count++;
    }
   
    assert1 (count == m, 15);

    // Check all the edges, and set the values up for CSA.
    for (int i = 0; i < m; i++) {
      assert1(igraph(i,0) >= 0 && igraph(i,0) < n, 16);
      assert1(igraph(i,1) >= 0 && igraph(i,1) < n, 17);
      igraph(i,0) += 1;
      igraph(i,1) += 1+n;
    }

    mexPrintf("Computing perfect matching.\n");
    // Solve the assignment problem.
    CSA csa(2*n,m,igraph.data());
    assert1(csa.edges()==n, 18);

    mexPrintf("Extracting solution from matching.\n");
    Array2D<int> ograph (n,3);
    for (int i = 0; i < n; i++) {
      int a,b,c;
      csa.edge(i,a,b,c);
      ograph(i,0)=a-1; ograph(i,1)=b-1-n; ograph(i,2)=c;
    }

    // Check the solution.
    // Count the number of high-cost edges from the perfect match
    // overlay that were used in the match.
    int overlayCount = 0;
    {
      std::pair<Label,Label> lp;
      std::pair<std::pair<Label,Label>, unsigned int> lpm;
      for (int a = 0; a < n; a++) {
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        const int c = ograph(a,2);
        assert1 (i >= 0 && i < n, 19);
        assert1 (j >= 0 && j < n, 20);
        // assert1 (c >= 0, 21);
        // edge from high-cost perfect match overlay
        if (c == ow * multiplier)
          overlayCount++;
        // skip outlier edges
        if (i >= n_label_0) { continue; }
        if (j >= n_label_1) { continue; }
        lp.first = node_id_2_label_0[i];
        lp.second = node_id_2_label_1[j];
        lpm.first = lp;
        lpm.second = c;
        matched_pair.push_back(lpm);
//        mexPrintf("%d %d %d\n", i, j, c);
      }
    }
    
    // Print a warning if any of the edges from the perfect match overlay
    // were used.  This should happen rarely.  If it happens frequently,
    // then the outlier connectivity should be increased.
    if (overlayCount > 0) {
      mexPrintf("WARNING: The match includes outlier(s) from the perfect match overlay.\n");
    }
    mexPrintf("done.\n");
  }

  plhs[0] = mxCreateNumericMatrix(matched_pair.size(), 3, mxUINT32_CLASS,
                                 mxREAL);
  {
    std::vector<std::pair<std::pair<Label, Label>, unsigned int> >::
      iterator mp_it;
    std::pair<Label,Label> lp;
    unsigned int i;
    int s = matched_pair.size();
    unsigned int * mp_out = (unsigned int *) mxGetPr(plhs[0]);
    for(mp_it=matched_pair.begin(), i=0; mp_it!=matched_pair.end();
        mp_it++, i++){
      lp = (*mp_it).first;
      
      mp_out[i] = lp.first;
      mp_out[i+s] = lp.second;
      mp_out[i+2*s] = (*mp_it).second;
    }
  }

  if(nlhs>=2){
    Label_Hash matched_0, matched_1;
    {
      std::vector<std::pair<std::pair<Label, Label>, unsigned int> >::
        iterator mp_it;
      std::pair<Label,Label> lp;
      Label l0, l1;
      double overlap_norm;
      for(mp_it=matched_pair.begin(); mp_it!=matched_pair.end();
          mp_it++){
        lp = (*mp_it).first;
        l0 = lp.first;
        l1 = lp.second;
        overlap_norm = (max_overlap - (double)(*mp_it).second)/
          ((double)std::min(n_pixel_label_0[l0], n_pixel_label_1[l1]));
//        mexPrintf("%d %d %d %d %d %f\n", l0, l1,
//                  n_pixel_label_0[l0], n_pixel_label_1[l1], (*mp_it).second,
//                  overlap_norm);
        if(overlap_norm < min_overlap_norm_threshold)
          continue;
        matched_0[l0] = l1;
        matched_1[l1] = l0;
      }
    }
    
    int ndim = 3;
    int dimensions[3] = {height, width, depth};
    plhs[1] = mxCreateNumericArray(ndim, dimensions, mxUINT8_CLASS, mxREAL);
    {
      unsigned char * d_out_0 = (unsigned char *) mxGetPr(plhs[1]);
      int i;
      unsigned char flag;
      for(i=0; i<n_pixel; i++){
        if(stack_0[i]<=0 || stack_1[i]<=0)
          continue;
        flag = 0;
        if(matched_0[stack_0[i]]==0){
          flag |= 0x1;
        }
        if(matched_1[stack_1[i]]==0){
          flag |= 0x4;
        }
        if(matched_1[stack_1[i]]!=0 && matched_0[stack_0[i]]!=0 &&
           matched_0[stack_0[i]]!=stack_1[i]){
          flag |= 0x2;
        }
        d_out_0[i] = flag;
      }
    }
  }
}
