// Some commonly used hash-functions
// 
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	09022009	init. code
//

#ifndef __HASH_FUNCTIONS__
#define __HASH_FUNCTIONS__

#include <ext/hash_map>

namespace std{
  using namespace __gnu_cxx;
}

struct eqi{ bool operator()(const int s1, const int s2) const
    {return s1==s2;}};
typedef  std::hash_map<int, int,
                       std::hash<int>, eqi> Hash_Int32_Int32;

struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int,
                       std::hash<unsigned int>, equi> Hash_UInt32_UInt32;

struct equli{ bool operator()(const unsigned long int s1, const unsigned long int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned long int, unsigned long int,
                       std::hash<unsigned long int>, equli> Hash_UInt64_UInt64;

typedef  std::hash_map<unsigned long int, unsigned int,
                       std::hash<unsigned long int>, equli> Hash_UInt64_UInt32;

struct eqd{ bool operator()(const double s1, const double s2) const
    {return s1==s2;}};
typedef  std::hash_map<double, double,
                       std::hash<double>, equi> Hash_Double_Double;

#endif
