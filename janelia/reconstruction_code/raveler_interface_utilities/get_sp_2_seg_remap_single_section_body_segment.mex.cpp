// Proofreaders may group segments together using the 3D link feature
// even though the body may be occuring in just one section and the
// proofreaders should have ideally used the add-segment tool. In such
// cases, it is not possible to import the proofread edits into the
// linkage graph (the body is in only one section). Instead, we must
// change the superpixel to segment map.
//
// Given the sp_2_seg and seg_2_body datastructures, return a
// remapped version of sp_2_seg.
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus,
// HHMI.
//

#include <mex.h>
#include <math.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <ext/hash_map>

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
typedef hash_map<unsigned long int, unsigned long int,
                 hash<unsigned long int>, equli> Hash_UInt64_UInt64;


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("input params:\n");
    mexPrintf("\t1.	superpixel-to-segment-map Px3 uint32\n");
    mexPrintf("\t2.	segment-to-body-map Qx2 uint32\n");
    mexPrintf("\tOutput:");
    mexPrintf("\t1. superpixel-to-segment-map Px3 uint32\n");
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

	const mxArray * sp_2_seg_mx = prhs[0];
  const mxArray * seg_2_body_mx = prhs[1];
  
  const int * size_1 = mxGetDimensions(sp_2_seg_mx);
  int n_sp_2_seg = std::max(size_1[0], size_1[1]);
  const int * size_2 = mxGetDimensions(seg_2_body_mx);
  int n_seg_2_body = std::max(size_2[0], size_2[1]);

  unsigned int * sp_2_seg = (unsigned int *) mxGetPr(sp_2_seg_mx);
  unsigned int * seg_2_body = (unsigned int *) mxGetPr(seg_2_body_mx);

  // create a segment-to-body hash
  Hash_UInt32_UInt32 seg_to_body_map;
  {
    int i;
    for(i=0; i<n_seg_2_body; i++)
      seg_to_body_map[seg_2_body[i]] = seg_2_body[i+n_seg_2_body];
  }

  // get the minimum and maximum section occurence for each body
  Hash_UInt32_UInt32 body_min_z, body_max_z;
  unsigned int z_max;
  unsigned int z_min;
  {
    // get the minimum and maximum section id
    int i;
    z_max = 0;
    {
      for(i=0; i<n_sp_2_seg; i++)
        z_max = std::max(z_max, sp_2_seg[i]);
    }
    z_min = z_max;
    {
      for(i=0; i<n_sp_2_seg; i++)
        z_min = std::min(z_min, sp_2_seg[i]);
    }

    // Set bounds for the minimum and maximum values
    {
      for(i=0; i<n_seg_2_body; i++){
        body_min_z[seg_2_body[i+n_seg_2_body]] = z_max;
        body_max_z[seg_2_body[i+n_seg_2_body]] = z_min;
      }
    }

    // now get the minimum and maximum z's
    unsigned int body_id;
    for(i=0; i<n_sp_2_seg; i++){
      body_id = seg_to_body_map[sp_2_seg[i+n_sp_2_seg*2]];
      body_min_z[body_id] = std::min(body_min_z[body_id], sp_2_seg[i]);
      body_max_z[body_id] = std::max(body_max_z[body_id], sp_2_seg[i]);
    }
  }

  // generate a remapping of segments. By default a segment gets
  // mapped to itself. But if the segment is getting mapped to a body
  // occuring in only one section then all the segments of such a body
  // get mapped to one segment. This will be used to remap the
  // superpixel-to-segment mapping.
  Hash_UInt32_UInt32 is_remapped_body;
  {
    Hash_UInt32_UInt32::iterator l_it;
    for(l_it=seg_to_body_map.begin(); l_it!=seg_to_body_map.end(); l_it++){
      if(body_min_z[(*l_it).second]==body_max_z[(*l_it).second]){
        // single section body
        if(is_remapped_body[(*l_it).second]==0){
          // First time seeing this body. Set the single segment to
          // which all other segments get mapped to
          is_remapped_body[(*l_it).second] = (*l_it).first;
        }
      }
    }
  }

  // output the remapped superpixel-to-segment-map
  plhs[0] = mxCreateNumericMatrix(n_sp_2_seg, 3, mxUINT32_CLASS, mxREAL);
  {
    unsigned int i;
    unsigned int * sp_2_seg_new = (unsigned int *) mxGetPr(plhs[0]);
    unsigned int seg_id;
    for(i=0; i<n_sp_2_seg; i++){
      sp_2_seg_new[i] = sp_2_seg[i];
      sp_2_seg_new[i+n_sp_2_seg] = sp_2_seg[i+n_sp_2_seg];
      seg_id = sp_2_seg[i+2*n_sp_2_seg];
      if(is_remapped_body[seg_to_body_map[seg_id]]==0)
        sp_2_seg_new[i+2*n_sp_2_seg] = seg_id;
      else
        sp_2_seg_new[i+2*n_sp_2_seg] = is_remapped_body[seg_to_body_map[seg_id]];
    }
  }
  
  return;
}
