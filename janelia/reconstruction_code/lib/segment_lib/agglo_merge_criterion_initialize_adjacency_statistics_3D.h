// macro routine
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI.
//

namespace iput
{
  namespace seg
  {
    bool PARENT_CLASS__::initialize_adjacency_statistics_3D__(
      const Array<Boundary_Value> * boundary_map,
      const Array<Label> * segment_map){ 
#ifdef __DEBUG__
      std::cout << "START: initialize_adjacency_statistics_3D__\n";
#endif
      int x, y, z;
      int width = boundary_map->dimensions[0];
      int height = boundary_map->dimensions[1];
      int depth = boundary_map->dimensions[2];
#ifdef __DEBUG__
      std::cout << "boundary_map: width: " << width << " height: " 
		<< height << " depth: " << depth << '\n';
      std::cout << "segment_map: width: " << segment_map->dimensions[0]
		<< " height: " << segment_map->dimensions[1]
		<< " depth: " << segment_map->dimensions[2] << '\n';
#endif
      int n_pixel_in_plane = width*height;
      Label * s_ptr;
      Label s0, s1, s2;
      unsigned char * b_ptr;

      for(z=0, s_ptr=segment_map->buffer, b_ptr=boundary_map->buffer;
	  z<depth; z++){
	for(y=0; y<height; y++){
	  for(x=0; x<width; x++, s_ptr++, b_ptr++){
	    s0 = *s_ptr;
	    if(s0==0){
	      if(y>0 && y<height-1){
		s1 = *(s_ptr-width);
		s2 = *(s_ptr+width);
		{COLLECT_BOUNDARY_STATISTICS__;};
	      }
	      if(x>0 && x<width-1){
		s1 = *(s_ptr-1);
		s2 = *(s_ptr+1);
		{COLLECT_BOUNDARY_STATISTICS__;};
	      }
	      if(z>0 && z<depth-1){
		s1 = *(s_ptr-n_pixel_in_plane);
		s2 = *(s_ptr+n_pixel_in_plane);
		{COLLECT_BOUNDARY_STATISTICS__;};
	      }
	    }
	    else{
	      {COLLECT_INTERNAL_STATISTICS__;};
	    }
	  }
	}
      }

      // prepare the adjacency list
      {
        Hash_Label_Pair_UInt32::iterator hp;
        for(hp=n_boundary_values.begin(); hp!=n_boundary_values.end(); hp++){
          Label_Pair_Key lpk = (*hp).first;
          adj_list.push_back(LPK_2_LP(lpk));
          Label s1, s2;
          boost::tie(s1,s2) = LPK_2_LP(lpk);
#ifdef __DEBUG__
          std::cout << "s1: " << s1 << " s2: " << s2 << '\n';
#endif
        }
      }

#ifdef __DEBUG__
      std::cout << "STOP: initialize_adjacency_statistics_3D__\n";
#endif
      return true;
    }
  }
}  
