// apply sequence of filters

#ifndef __IPUJ_FILTER__
#define __IPUJ_FILTER__

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <array.h>

namespace iput
{
  template <typename T>
    int filter(const Array<T> & source, Array<T> & target,
	       const std::string & filter_param, int is_verbose=0){
    if(is_verbose>0)
      std::cout << "START: filter\n";
    Array<T> temp_source = source;

    std::istringstream iss (filter_param, std::istringstream::in);
    while(iss.good()){
      target = source;

      std::string fn;
      iss >> fn;
      if(is_verbose>1)
	std::cout << "fn: " << fn << '\n';
      
      if(fn.compare("neg")==0){
	T p;
	{
	  double tp;
	  iss >> tp;
	  p = tp;
	  if(is_verbose>1)
	    std::cout << "p: " << tp << '\n';
	}
	for(int i=0; i<source.n_element; i++)
	  target.buffer[i] = p - temp_source.buffer[i];
      }

      if(fn.compare("LT")==0){
	T p;
	{
	  double tp;
	  iss >> tp;
	  p = tp;
	  if(is_verbose>1)
	    std::cout << "p: " << tp << '\n';
	}
	for(int i=0; i<source.n_element; i++)
	  target.buffer[i] = temp_source.buffer[i]<p?p:temp_source.buffer[i];
      }

      if(fn.compare("UT")==0){
	T p;
	{
	  double tp;
	  iss >> tp;
	  p = tp;
	  if(is_verbose>1)
	    std::cout << "p: " << tp << '\n';
	}
	for(int i=0; i<source.n_element; i++)
	  target.buffer[i] = temp_source.buffer[i]>p?p:temp_source.buffer[i];
      }

      temp_source = target;
    }
    if(is_verbose>0)
      std::cout << "STOP: filter\n";
    return 0;
  }
}

#endif
