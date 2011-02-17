// Given a linkage graph generate segment-to-body-map. Negative
// segment ids are for dummy segments not to be included in the
// segment-to-body-map.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <ext/hash_map>

#include <merge_sets_h.h>

using namespace __gnu_cxx;

struct equi{
  bool operator()(const unsigned int s1, const unsigned int s2) const
    {
      return s1==s2;
    }
};
typedef hash_map<unsigned int, unsigned int,
                 hash<unsigned int>,
                      equi> Hash_UInt32_UInt32;

struct equli{
  bool operator()(const unsigned long int s1, const unsigned long int s2) const
    {
      return s1==s2;
    }
};

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	linkage graph\n");
    mexPrintf("\t2. linkage threshold\n");
    mexPrintf("\tOutput:");
    mexPrintf("\t1. segment-to-body map Rx1 uint32\n");
    return;
  }
  if(nrhs!=2)
  {
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }

  if(nlhs>1)
  {
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  mexPrintf("START: get_sec_seg_2_body_map_from_links_3D_with_dummy\n");
	const mxArray * linkage_graph_mx = prhs[0];
  const mxArray * linkage_threshold_mx = prhs[1];
  
  const int * size_1 = mxGetDimensions(linkage_graph_mx);
  int n_plane = std::max(size_1[0], size_1[1]);
  double linkage_threshold = * mxGetPr(linkage_threshold_mx);
  mexPrintf("linkage_threshold: %g\n", linkage_threshold);
  
  const int section_bit_offset = 32;
  const unsigned long int segment_mask = 0x00000000ffffffff;
  
  // Build Merge_Sets_H datastructure for [section_id, segment_id]
  // tuples.
  Merge_Sets_H<unsigned long int, hash<unsigned long int>, equli> M(NULL);
  {
    double * section_pair_links;
    int i;
    unsigned long int l1, l2, l;
    for(i=0; i<n_plane; i++){
      mexPrintf("plane: %d\n", i);
      section_pair_links = mxGetPr(mxGetCell(linkage_graph_mx, i));
      const int * size_l = mxGetDimensions(mxGetCell(linkage_graph_mx, i));
      int n_pair = size_l[0];
      mexPrintf("n_pair: %d\n", n_pair);
      int j;
      for(j=0; j<n_pair; j++){
        l1 = i;
        l1 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j];
        l1 |= l;
        M.add_new_set_inc(l1);

        l2 = i+1;
        l2 <<= section_bit_offset;
        l = (unsigned int) (int) section_pair_links[j+n_pair];
        l2 |= l;
        M.add_new_set_inc(l2);

        if(section_pair_links[j]!=0 &&
           section_pair_links[j+n_pair]!=0 &&
           section_pair_links[j+2*n_pair]>linkage_threshold){
          M.merge(l1, l2);
        }
        
        // if dummy node then add an extra node in the next section in
        // case the body disappears for two sections.
        if(section_pair_links[j]!=0 &&
           section_pair_links[j+n_pair]<0 &&
           section_pair_links[j+2*n_pair]>linkage_threshold){
          l2 = i+2;
          l2 <<= section_bit_offset;
          l = (unsigned int) (int) section_pair_links[j+n_pair];
          l2 |= l;
          M.add_new_set_inc(l2);

          M.merge(l1, l2);
        }
      }
    }
  }

  // Construct segment-to-body-map
  {
    // get the number of valid segments and bodies.
    unsigned int n_segment_real = 0;
    Hash_UInt32_UInt32 body_id_hash;
    unsigned int a, n, p;
    unsigned long int l, l1;
    int s;
    unsigned int max_p=0;
    Merge_Sets_H<unsigned long int, hash<unsigned long int>, equli>::
      Adam_Hash_Iterator it;
    for(it=M.adam.begin(); it!=M.adam.end(); it++){
      l = (*it).first;
      l1 = l & segment_mask;
      s = (int) l1;
      if(s<=0)
        continue;
      n_segment_real++;

      l1 = l >> section_bit_offset;
      p = (unsigned int) l1;
      max_p = std::max(max_p, p);
      
      a = M.get_adam_id(l);
      if(body_id_hash[a]==0){
        n = body_id_hash.size();
        body_id_hash[a] = n;
      }
    }
    mexPrintf("n_segment_real: %d\n", n_segment_real);
    mexPrintf("max plane: %u\n", max_p);
    
    plhs[0] = mxCreateNumericMatrix(n_segment_real, 3, mxUINT32_CLASS, mxREAL);
    unsigned int * sec_seg_2_body = (unsigned int *) mxGetPr(plhs[0]);
    for(it=M.adam.begin(), n=0; it!=M.adam.end(); it++){
      l = (*it).first;
      l1 = l & segment_mask;
      s = (int) l1;
      if(s<=0)
        continue;
      l1 = l >> section_bit_offset;
      p = (unsigned int) l1;
      sec_seg_2_body[n] = p; // plane or layer id
      sec_seg_2_body[n+n_segment_real] = s; // segment id
      sec_seg_2_body[n+2*n_segment_real] = body_id_hash[M.get_adam_id(l)];
      
      n++;
    }
  }
  mexPrintf("STOP: get_sec_seg_2_body_map_from_links_3D_with_dummy\n");
  return;
}
