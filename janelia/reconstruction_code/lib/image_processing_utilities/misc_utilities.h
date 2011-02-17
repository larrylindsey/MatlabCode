// miscellaneous utitlity functions
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//

#ifndef __IPUJ_MISC_UTILITIES__
#define __IPUJ_MISC_UTILITIES__

#include <vector>
#include <string>

namespace iput
{
  std::vector<std::string> split_param_string(const std::string & param, 
					      char split_char,
					      char escape_char);

  std::vector<int> find_char_in_string(const std::string & s, 
				       char query_char,
				       char escape_char);
}

#endif // __IPUJ_MISC_UTILITIES__
