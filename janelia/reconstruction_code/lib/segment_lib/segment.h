// base class for segmentation
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

#ifndef __IPUJ_SEGMENT__
#define __IPUJ_SEGMENT__

#include <ext/hash_map>

#include <array.h>

namespace std{
  using namespace __gnu_cxx;
}

namespace iput
{
namespace seg
{
  typedef unsigned char Boundary_Value; // has to be integral
  typedef unsigned int Label;
  struct eqlp{ bool operator()(const Label s1, const Label s2) const
      {return s1==s2;}};
  // Hash maps of labels
  typedef std::hash_map<Label, unsigned int, std::hash<Label>,
                        eqlp> Hash_Label_UInt32;


  typedef std::pair<Label,Label> Label_Pair;
  typedef unsigned long int Label_Pair_Key;
#define LP_2_LPK(lp) (lp.first<lp.second?((((unsigned long int)lp.first)<< 32) | ((unsigned long int) lp.second)):((((unsigned long int)lp.second)<< 32) | ((unsigned long int) lp.first)))
#define LPK_2_LP(lpk) (iput::seg::Label_Pair((iput::seg::Label)(lpk>>32), (iput::seg::Label)(lpk & 0x00000000ffffffff)))
  struct eqlpk{ bool operator()(const Label_Pair_Key s1,
                               const Label_Pair_Key s2) const
      {return s1==s2;}};
  // Hash maps of label pairs
  typedef std::hash_map<Label_Pair_Key, unsigned int, std::hash<Label_Pair_Key>
                        , eqlpk> Hash_Label_Pair_UInt32;
  
  
  class Segment {
  public:
    // functions to be declared by derived classes
    virtual bool initialize(const Array<Boundary_Value> * boundary_map, 
			    const Array<Label> * initial_segment_map) = 0;
    virtual bool compute_segmentation(void) = 0;

    int is_verbose;
    
    // functions to be called by derived classes
    inline Label_Pair_Key lp_2_lpk(Label_Pair lp){
      return LP_2_LPK(lp);
    }
  };
}
}
#endif //__IPU_SEGMENT__

