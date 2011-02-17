// split a string of sequence of parameters into individual parameters


#include <iostream>
#include <vector>
#include <string>

#include <misc_utilities.h>

namespace iput
{
  std::vector<std::string> split_param_string(const std::string & param, 
				      char split_char,
				      char escape_char){
    std::vector<std::string> sub_params;

    std::vector<int> splits = 
      find_char_in_string(param, split_char, escape_char);

    int prev_split = -1;
    for(int i=0; i<splits.size(); i++){
      if(splits[i]>prev_split+1){
	sub_params.push_back(param.substr(prev_split+1, 
					     splits[i]-prev_split-1));
      }
      else
	sub_params.push_back(std::string(""));

      prev_split = splits[i];
    }

    sub_params.push_back(param.substr(prev_split+1, 
				  param.length() - prev_split));

    return sub_params;
  }

  std::vector<int> find_char_in_string(const std::string & s, 
				       char query_char,
				       char escape_char){
    std::vector<int> l;
    int i;
    int state = 1;
    for(i=0; i<s.length(); i++){
      if(state==1){
	if(s[i]=='\0')
	  break;
	if(s[i]==query_char)
	  l.push_back(i);
	if(s[i]==escape_char)
	  state=2;
	continue;
      }
      if(state==2){
	state=1;
	continue;
      }
    }
    return l;
  }
}

