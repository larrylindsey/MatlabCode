// testbed for agglomerative segmentation on region adjacency graph
// using mean boundary ladder criterion
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include <image_lib.h>
#include <utilities.h>

#include <array.h>
#include <filter.h>
#include <hash_functions.h>
#include <file_input_output.h>
#include <segment_agglo_RAG.h>
#include <agglo_merge_criterion.h>
#include <agglo_merge_criterion_mean_boundary.h>
#include <agglo_merge_criterion_mean_boundary_ladder.h>

int main(int argc, char * argv[]){

  static char * Spec [] = { "[!V] [-v <int>] [-{a|min_area_threshold} <int>]",
			    " [-{f|filter} <string>]",
			    " [-{SW|shrink_wrap}] ",
			    " <boundary_map:string>", 
			    " <initial_segment_map:string>",
			    " [-{s|intermediate_segment_map} <string>]",
			    " <output_segment_map:string>",
			    " <f_thresholds:string>",
			    NULL  };

  Process_Arguments(argc, argv, Spec, 0);

  int is_verbose = 0;
  if(Is_Arg_Matched("-V"))
    is_verbose = 255;
  if(Is_Arg_Matched("-v @"))
    is_verbose = Get_Int_Arg("-v @");
  
  if(is_verbose>0)
    std::cout << "START: " << argv[0] << '\n';

  if(is_verbose>1)
    std::cout << "Reading in boundary map\n";
  Stack * stack = Read_Stack(Get_String_Arg("boundary_map"));
  iput::Array<unsigned char> boundary_map;
  boundary_map.n_dimension = 3;
  boundary_map.dimensions.push_back(stack->width);
  boundary_map.dimensions.push_back(stack->height);
  boundary_map.dimensions.push_back(stack->depth);
  boundary_map.n_element = stack->width*stack->height*stack->depth;
  boundary_map.buffer = (unsigned char *) stack->array;
  if(boundary_map.buffer==NULL){
    printf("Could not read boundary map.\n");
    return -1;
  }

  if(Is_Arg_Matched("-f @")){
    char * filter_param = Get_String_Arg("-f @");
    if(is_verbose>1)
      std::cout << "Applying filter " << filter_param << '\n';
    iput::Array<unsigned char> b;
    iput::filter(boundary_map, b, filter_param, is_verbose-1);
    boundary_map = b;
  }

  if(is_verbose>1)
    std::cout << "Reading in initial segment map\n";
  iput::Array<unsigned int> initial_segment_map =
    iput::Array<unsigned int>::
    fread_raw_array(Get_String_Arg("initial_segment_map"));
  if(initial_segment_map.buffer==NULL){
    printf("Could not read initial segment map.\n");
    return -1;
  }

  iput::Array<unsigned int> interm_segment_map;
  if(Is_Arg_Matched("-s @")){
    if(is_verbose>1)
      std::cout << "Reading in intermediate segment map\n";
    interm_segment_map = 
      iput::Array<unsigned int>::fread_raw_array(Get_String_Arg("-s @"));
    if(interm_segment_map.buffer==NULL){
      printf("Could not read intermediate segment mapping.\n");
      return -1;
    }

    Hash_UInt32_UInt32 im;
    int i;
    for(i=0; i<interm_segment_map.dimensions[1]; i++)
      im[interm_segment_map.buffer[i*2]] = 
	interm_segment_map.buffer[i*2+1];
    int n_voxel = stack->width*stack->height*stack->depth;
    for(i=0; i<n_voxel; i++)
      initial_segment_map.buffer[i] = 
	im[initial_segment_map.buffer[i]];
  }

  bool is_seeded_seg = false;
  if(Is_Arg_Matched("-SW"))
    is_seeded_seg = true;

  if(is_verbose>1)
    std::cout << "Setting up merge criterion\n";
  iput::seg::Segment_Agglo_RAG S;
  S.is_verbose = is_verbose - 1;
  iput::seg::Agglo_Merge_Criterion * MC;
  if(Is_Arg_Matched("-a @")){
    if(is_verbose>1)
      std::cout << "agglo mean boundary ladder\n";
    iput::seg::Agglo_Merge_Criterion_Mean_Boundary_Ladder * M =
      new iput::seg::Agglo_Merge_Criterion_Mean_Boundary_Ladder();
    M->area_threshold = Get_Int_Arg("-a @");
    if(is_verbose>1)
      std::cout << "M->area_threshold: " << M->area_threshold << '\n';
    MC = (iput::seg::Agglo_Merge_Criterion *) M;
    S.merge_criterion = MC;
    S.is_enabled_break_at_last_f_threshold = false;
  }
  else{
    if(is_verbose>1)
      std::cout << "agglo mean boundary\n";
    MC = (iput::seg::Agglo_Merge_Criterion *) 
      new iput::seg::Agglo_Merge_Criterion_Mean_Boundary();
    S.merge_criterion = MC;
    S.is_enabled_break_at_last_f_threshold = true;
  }
  MC->is_verbose = is_verbose - 1;

  if(is_seeded_seg){
    S.is_enabled_break_at_last_f_threshold = false;
    if(is_verbose>1)
      std::cout << "Preparing seeds for shrink wrap\n";
    Hash_UInt32_UInt32 is_seeded;
    { // top and bottom faces
      int offset = (stack->depth-1)*stack->width*stack->height;
      for(int i=0; i<stack->width*stack->height; i++){
	is_seeded[initial_segment_map.buffer[i]]=1;
	is_seeded[initial_segment_map.buffer[i+offset]]=1;
      }
    }
    { // east and west faces
      int np = stack->width*stack->height;
      for(int d=0; d<stack->depth; d++)
	for(int y=0; y<stack->height; y++){
	  is_seeded[initial_segment_map.buffer[d*np+y*stack->width]]=1;
	  is_seeded[initial_segment_map.buffer[d*np+y*stack->width+
					       stack->width-1]]=1;
	}
    }
    { // north and south faces
      int np = stack->width*stack->height;
      for(int d=0; d<stack->depth; d++)
	for(int x=0; x<stack->width; x++){
	  is_seeded[initial_segment_map.buffer[d*np+x]]=1;
	  is_seeded[initial_segment_map.buffer[d*np+x+np-stack->width]]=1;
	}
    }
    if(is_verbose>16)
      std::cout << "Seed segment labels:\n";
    Hash_UInt32_UInt32::iterator i;
    for(i=is_seeded.begin(); i!=is_seeded.end(); i++)
      if((*i).second!=0){
	if(is_verbose>16)
	  std::cout << (*i).first << ' ';
	S.seeds.push_back((*i).first);
      }
    if(is_verbose>16)
      std::cout << '\n';
  }

  if(is_verbose>1)
    std::cout << "Initializing segmentation\n";
  if(!S.initialize(& boundary_map, & initial_segment_map))
    return -1;
  else
    if(is_verbose>1)
      std::cout << "Segmentation initialization done.\n";

  {
    int i;
    char * f_s = Get_String_Arg("f_thresholds");
    std::stringstream f_ss(f_s, std::ios::in);
    while(f_ss.good()){
      double f;
      f_ss >> f;
      S.f_thresholds.push_back(f);
    }
  }

  if(is_verbose>1)
    std::cout << "Computing segmentation\n";
  {
    bool result;
    if(is_seeded_seg)
      result = S.compute_segmentation_seeded();
    else
      result = S.compute_segmentation();
    
    if(!result)
      return -1;
    else
      if(is_verbose>1)
	std::cout << "Segmentation done.\n";
  }

  {
    char * filename = (char *) new 
      char[strlen(Get_String_Arg("output_segment_map"))+100];
    for(int j=0; j<S.f_thresholds.size(); j++){
      if(is_verbose>2)
	std::cout << "Saving for f_threshold: " << S.f_thresholds[j] << '\n';

      if(S.segment_mappings[j].n_dimension==-1)
	continue;

      {
	char sf[100];
	sprintf(sf, "%g", S.f_thresholds[j]);
	sprintf(filename, Get_String_Arg("output_segment_map"), sf);
      }
      if(is_verbose>1)
	std::cout << "output filename: " << filename << '\n';

      if(Is_Arg_Matched("-s @")){ // map the initial segments to the new ones
	if(is_verbose>1)
	  std::cout << "mapping the initial segments to the new ones\n";
	Hash_UInt32_UInt32 nm;
	int i;
	for(i=0; i<S.segment_mappings[j].dimensions[1]; i++)
	  nm[S.segment_mappings[j].buffer[i*2]] = 
	    S.segment_mappings[j].buffer[i*2+1];
	iput::Array<unsigned int> nma;
	nma = interm_segment_map;
	for(i=0; i<nma.dimensions[1]; i++)
	  nma.buffer[i*2+1] = nm[nma.buffer[i*2+1]];
	nma.fwrite_raw_array(filename);
      }
      else{
	S.segment_mappings[j].fwrite_raw_array(filename);
      }
      
    }
    
    if(!S.is_enabled_break_at_last_f_threshold){
      {
	char sf[100];
	sprintf(sf, "L%g", S.f_thresholds.back());
	sprintf(filename, Get_String_Arg("output_segment_map"), sf);
      }
      if(is_verbose>1)
	std::cout << "output filename: " << filename << '\n';
      if(Is_Arg_Matched("-s @")){ // map the initial segments to the new ones
	if(is_verbose>1)
	  std::cout << "mapping the initial segments to the new ones\n";
	Hash_UInt32_UInt32 nm;
	int i;
	int j = S.f_thresholds.size();
	for(i=0; i<S.segment_mappings[j].dimensions[1]; i++)
	  nm[S.segment_mappings[j].buffer[i*2]] = 
	    S.segment_mappings[j].buffer[i*2+1];
	iput::Array<unsigned int> nma;
	nma = interm_segment_map;
	for(i=0; i<nma.dimensions[1]; i++)
	  nma.buffer[i*2+1] = nm[nma.buffer[i*2+1]];
	nma.fwrite_raw_array(filename);
      }
      else{
	S.segment_mappings[S.f_thresholds.size()].fwrite_raw_array(filename);
      }
    }

    delete [] filename;
  }

  delete MC;
  Free_Stack(stack);

  if(is_verbose>0)
    std::cout << "STOP: " << argv[0] << '\n';
  return 0;
}

