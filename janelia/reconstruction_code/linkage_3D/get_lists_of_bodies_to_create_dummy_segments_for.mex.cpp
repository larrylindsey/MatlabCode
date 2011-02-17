// Given the superpixel-to-segment mapping [z, sp_id, seg_id] and
// segment-to-body mapping [seg_id, body_id], find discontinuities in
// the bodies. Such discontinuities may occur if a body "disappears"
// entirely from a section, e.g., due to a fold.  Generates a cell
// array, with one cell for each section. Each cell contains the list
// of bodies "disappearing" in that section.
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
    mexPrintf("\t1. cell array Zx1 of uint32 vectors\n");
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

  mexPrintf("z_min: %u\n", z_min);
  mexPrintf("z_max: %u\n", z_max);

  // create a hash for [body_id, section_id] pairs
  Hash_UInt64_UInt64 body_section_pair;
  {
    unsigned long int l1, l;
    int i;
    for(i=0; i<n_sp_2_seg; i++){
      l1 = seg_to_body_map[sp_2_seg[i+n_sp_2_seg*2]];
      l1 <<= 32;
      l = sp_2_seg[i];
      l |= l1;
      body_section_pair[l] = 1;
    }
  }
  
  //
  // create a list for each section of the disappearing body ids
  //
  std::vector<std::vector<unsigned int> > lists_disappearing_body;
  // initialize
  {
    int i;
    std::vector<unsigned int> empty_list;
    for(i=0; i<=z_max-z_min; i++)
      lists_disappearing_body.push_back(empty_list);
  }

  // loop through bodies and find disappearing cases
  {
    Hash_UInt32_UInt32::iterator it;
    unsigned int body_id;
    unsigned int z1, z2;
    unsigned int z;
    unsigned long int l1, l;
    for(it=body_min_z.begin(); it!=body_min_z.end(); it++){
      body_id = (*it).first;
      z1 = (*it).second;
      z2 = body_max_z[body_id];
      for(z=z1; z<=z2; z++){
        l1 = body_id;
        l1 <<= 32;
        l = z;
        l |= l1;
        if(body_section_pair[l]==0){
          lists_disappearing_body[z-z_min].push_back(body_id);
        }
      }
    }
  }

  // output to cell array
  plhs[0] = mxCreateCellMatrix(z_max, 1);
  {
    std::vector<std::vector<unsigned int> >::iterator it;
    std::vector<unsigned int> list;
    std::vector<unsigned int>::iterator it_l;
    mxArray * temp_list_mx;
    unsigned int * temp_list;
    unsigned int z, i;
    for(it=lists_disappearing_body.begin(), z=z_min;
        it!=lists_disappearing_body.end(); it++, z++){
      list = (*it);
      
      temp_list_mx = mxCreateNumericMatrix(list.size(), 1, mxUINT32_CLASS,
                                           mxREAL);
      
      temp_list = (unsigned int *) mxGetPr(temp_list_mx);
      for(it_l=list.begin(), i=0; it_l!=list.end(); it_l++, i++)
        temp_list[i] = *it_l;
      
      mxSetCell(plhs[0], z-1, temp_list_mx);
    }
  }
  return;
}
