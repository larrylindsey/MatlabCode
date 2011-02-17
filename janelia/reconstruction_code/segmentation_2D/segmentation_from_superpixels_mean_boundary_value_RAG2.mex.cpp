// compute segmentation from a superpixel map (agglomerative clustering)
// Feature: between each segment, compute the mean boundary value and merge them if this value is below
// a specified threshold
// 
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	08042008	init. code
//

#include <mex.h>
#include <math.h>
#include <stdio.h>
#include <ext/hash_map>
#include <utility>                       // for std::pair
#include <algorithm>                     // for std::for_each
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include <merge_sets_h.h>

#define M_PI_2 1.57079632679489661923
#define M_PI 3.14159265358979323846

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

namespace std{
  using namespace __gnu_cxx;
}

// Hash function for unsigned int to unsigned int
struct equi{ bool operator()(const unsigned int s1, const unsigned int s2) const
    {return s1==s2;}};
typedef  std::hash_map<unsigned int, unsigned int, std::hash<unsigned int>, equi> Hash_UInt32;
typedef  std::hash_map<unsigned int, int, std::hash<unsigned int>, equi> Hash_UInt32_Int32;

// priority queue for unsigned int
struct ltui{ bool operator()(const unsigned int f1, const unsigned int f2){ return f1<f2; }};
typedef std::map<const unsigned int, unsigned int, ltui> Q_UInt32;

// Hash function from pair of labels to unsigned int
typedef unsigned long int Label_Pair;
struct equli{ bool operator()(const Label_Pair s1, const Label_Pair s2) const
    {return s1==s2;}};
typedef  std::hash_map<Label_Pair, unsigned int, std::hash<Label_Pair>, equli> Label_Pair_Hash;

Label_Pair label_pair_2_id(unsigned int s1, unsigned int s2){
  Label_Pair lp1;
  lp1 = s1;
  lp1 <<= 32;
  Label_Pair lp;
  lp = s2;
  lp |= lp1;
  return lp;
}
std::pair<unsigned int, unsigned int> id_2_label_pair(Label_Pair lp){
  std::pair<unsigned int, unsigned int> s;
  s.first = lp >> 32;
  s.second = lp & 0x00000000ffffffff;
  return s;
}

// priority queue for label pairs
struct ltuli{ bool operator()(const Label_Pair f1, const Label_Pair f2){ return f1<f2; }};
typedef std::map<const Label_Pair, unsigned int, ltuli> Q_Label_Pair_UInt32;

// Priority queue for merging
struct ltdb{ bool operator()(const double f1, const double f2){ return f1<f2; }};
typedef std::map<const double, Label_Pair, ltdb> Q_Double_Label_Pair;

typedef std::hash_map<Label_Pair, Q_Double_Label_Pair::iterator, std::hash<Label_Pair>, equli> Hash_Label_Pair_Q_It;


using namespace boost;

typedef struct {
  unsigned int neighbor_id;
  int next_ptr;
} Adjacency_List_Node;

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		if(nlhs==1){
			plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
			*mxGetPr(plhs[0]) = 1;
			return;
		}
		mexPrintf("Usage old_label_to_new_label = segmentation_from_superpixels_mean_boundary_value_RAG(superpixel_label_map, boundary, f_threshold_seq, boundary_length_seq);\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. superpixel label map MxN matrix (uint32)\n");
		mexPrintf("\t2. scalar field on which watershed was performed MxN matrix (uint8)\n");
		mexPrintf("\t3. f_thresholds for which segmentation is to be saved 1xR (double).\n");
		mexPrintf("\t4. minimum length of boundaries to be considered for merging 1xR (uint32).\n");    
		mexPrintf("output:\n");
		mexPrintf("\t1. mapping from input segment labels to merged segment labels 1xR cell array of 1xS uint32 (S=max. label id).\n");
		return;
	}
	if(nrhs!=4){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>1){
		mexErrMsgTxt("Wrong number of outputs\n");
		return;
	}

  mexPrintf("START: segmentation_from_superpixels_mean_boundary_value_RAG2\n");

#ifdef __DEBUG__
  FILE * fout_db = fopen("/groups/chklovskii/home/vitaladevunis/temp/debug.log", "wt");
  fprintf(fout_db, "START: segmentation_from_superpixels_mean_boundary_value_RAG2\n");
#endif //__DEBUG__
  
  int numDim = mxGetNumberOfDimensions(prhs[0]);
	if(numDim!=2){
		mexErrMsgTxt("Wrong no. of dimensions for arg. 1\n");
		return;
	}

	const int * sizeImage;
	sizeImage = mxGetDimensions(prhs[0]);
	int width = sizeImage[1], height = sizeImage[0];
	int n_pixel = height*width;
	unsigned int * superpixel_label_map = (unsigned int *) mxGetPr(prhs[0]);
	
	unsigned char * boundary_map = (unsigned char *) mxGetPr(prhs[1]);
	double * f_thresholds = mxGetPr(prhs[2]);
  int n_f_threshold = mxGetN(prhs[2]);
	unsigned int * boundary_length_thresholds = (unsigned int *) mxGetPr(prhs[3]);

  Label_Pair_Hash sum_boundary_value;
  Label_Pair_Hash n_boundary_pixel;
  // get sum and number of boundary values for each pair of adjacent segments
#ifdef __PROFILE__
  mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
  Label_Pair_Hash label_adjacency;
  unsigned int max_superpixel_id = 0;
  Hash_UInt32 n_neighbor;
  Hash_UInt32_Int32 adjacency_list_pointer;
  Adjacency_List_Node * adjacency_lists;
  Q_Double_Label_Pair merge_q;
  Hash_Label_Pair_Q_It merge_q_pos;
  Label_Pair_Hash is_valid_merge_q_pos;
  Merge_Sets_H<unsigned int, std::hash<unsigned int>, equi>
    superpixel_sets(NULL);
  {
    Q_Label_Pair_UInt32 label_pair_list;
    {
      int x, y;
      unsigned int * s;
      unsigned char * b;
      Label_Pair lp;
      for(x=1, s=superpixel_label_map+height, b=boundary_map+height;
          x<width-1;
          x++, s+=height, b+=height){
        if(*s==0){
          lp = label_pair_2_id(*(s-height), *(s+height));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
          lp = label_pair_2_id(*(s+height), *(s-height));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        else{
          max_superpixel_id = MAX(max_superpixel_id, *s);
          superpixel_sets.add_new_set_inc(*s);
        }
      }
      for(x=1, s=superpixel_label_map+height*2-1, b=boundary_map+height*2-1;
          x<width-1;
          x++, s+=height, b+=height){
        if(*s==0){
          lp = label_pair_2_id(*(s-height), *(s+height));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
          lp = label_pair_2_id(*(s+height), *(s-height));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        else{
          max_superpixel_id = MAX(max_superpixel_id, *s);
          superpixel_sets.add_new_set_inc(*s);
        }
      }
      for(y=1, s=superpixel_label_map+1, b=boundary_map+1;
          y<height-1;
          y++, s++, b++){
        if(*s==0){
          lp = label_pair_2_id(*(s-1), *(s+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
          lp = label_pair_2_id(*(s+1), *(s-1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        else{
          max_superpixel_id = MAX(max_superpixel_id, *s);
          superpixel_sets.add_new_set_inc(*s);
        }
      }
      for(y=1, s=superpixel_label_map+height*(width-1)+1, b=boundary_map+height*(width-1)+1;
          y<height-1;
          y++, s++, b++){
        if(*s==0){
          lp = label_pair_2_id(*(s-1), *(s+1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
          lp = label_pair_2_id(*(s+1), *(s-1));
          sum_boundary_value[lp] += *b;
          n_boundary_pixel[lp] ++;
        }
        else{
          max_superpixel_id = MAX(max_superpixel_id, *s);
          superpixel_sets.add_new_set_inc(*s);
        }
      }
      for(x=1, s=superpixel_label_map+height, b=boundary_map+height; // skip first column
          x<width-1;
          x++){
        s++; // skip first element in column
        b++;
        for(y=1;
            y<height-1;
            y++, s++, b++){
          if(*s==0){
            lp = label_pair_2_id(*(s-height), *(s+height));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
            lp = label_pair_2_id(*(s+height), *(s-height));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
          }
          else{
            max_superpixel_id = MAX(max_superpixel_id, *s);
            superpixel_sets.add_new_set_inc(*s);
          }
          if(*s==0){
            lp = label_pair_2_id(*(s-1), *(s+1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
            lp = label_pair_2_id(*(s+1), *(s-1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
          }
          if(*s==0){
            lp = label_pair_2_id(*(s-height-1), *(s+height+1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
            lp = label_pair_2_id(*(s+height+1), *(s-height-1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
          }
          if(*s==0){
            lp = label_pair_2_id(*(s-height+1), *(s+height-1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
            lp = label_pair_2_id(*(s+height-1), *(s-height+1));
            sum_boundary_value[lp] += *b;
            n_boundary_pixel[lp] ++;
          }
        }
        s++; // skip last element in column.
        b++;
      }

      // Build a sorted list of label pairs
      {
        Label_Pair lp;
        Label_Pair_Hash::iterator lp_it;
        for(lp_it=sum_boundary_value.begin(); lp_it!=sum_boundary_value.end();
            lp_it++){
          lp = (*lp_it).first;
          label_pair_list[lp] = 1;
          label_adjacency[lp] = 1;
        }
      }
#ifdef __DEBUG2__
      {
        fprintf(fout_db, "Label pairs:\n");
        Q_Label_Pair_UInt32::iterator lpq_it;
        unsigned int s0, s1;
        for(lpq_it=label_pair_list.begin(); lpq_it!=label_pair_list.end();
            lpq_it++){
          lp = (*lpq_it).first;
          tie(s0,s1) = id_2_label_pair(lp);
          fprintf(fout_db, "label0: %u, label1: %u\n", s0, s1);
        }
        fprintf(fout_db, "Boundary statistics 1:\n");
        Label_Pair_Hash::iterator lp_it;
        for(lp_it=sum_boundary_value.begin(); lp_it!=sum_boundary_value.end();
            lp_it++){
          lp = (*lp_it).first;
          tie(s0,s1) = id_2_label_pair(lp);
          fprintf(fout_db, 
            "label0: %u, label1: %u, sum_boundary_value: %u\n",
            s0, s1, (*lp_it).second);
        }
        fprintf(fout_db, "Boundary statistics 2:\n");
        for(lp_it=n_boundary_pixel.begin(); lp_it!=n_boundary_pixel.end();
            lp_it++){
          lp = (*lp_it).first;
          tie(s0,s1) = id_2_label_pair(lp);
          fprintf(fout_db, 
            "label0: %u, label1: %u, n_boundary_pixel: %u\n",
            s0, s1, (*lp_it).second);
        }
      }
#endif // __DEBUG2__
    }
#ifdef __PROFILE__
    mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif

    // Construct region adjacency graph with segments as nodes and edges between
    // spatially adjacent segments.
#ifdef __PROFILE__
    mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
    {
      adjacency_lists =
        (Adjacency_List_Node *) new Adjacency_List_Node[label_pair_list.size()+1];
      Q_Label_Pair_UInt32::iterator lpq_it;
      unsigned int s0, s1;
      unsigned int n, b;
      int node_id=1;
      int prev_node_id=-1;
      Label_Pair lp;
      double m;
      Q_Double_Label_Pair::iterator q_it;
      for(lpq_it=label_pair_list.begin(); lpq_it!=label_pair_list.end();
          lpq_it++){
        lp = (*lpq_it).first;
        tie(s0,s1) = id_2_label_pair(lp);
        if(s0==0 || s1==0 || s0==s1)
          continue;
        if(adjacency_list_pointer[s0]==0){
          adjacency_list_pointer[s0] = node_id;
          prev_node_id = -1;
        }
        if(prev_node_id!=-1)
          adjacency_lists[prev_node_id].next_ptr = node_id;
        
        n_neighbor[s0]++;
        
        adjacency_lists[node_id].next_ptr = -1;
        adjacency_lists[node_id].neighbor_id = s1;
        
        n = n_boundary_pixel[lp];
        b = sum_boundary_value[lp];
        m = (double)b/(double)n; // mean boundary value
        // insert in merge queue
        q_it = merge_q.find(m); // keep priority keys unique.
        while(q_it!=merge_q.end()){
          m += 0.00001*((double)(rand()%1000)); // epsilon increments
          q_it = merge_q.find(m);
        }
        {
          const std::pair<double, Label_Pair> v(m, lp);
          bool is_inserted;
          tie(q_it, is_inserted)  = merge_q.insert(v);
        }
        merge_q_pos[lp] = q_it;
        is_valid_merge_q_pos[lp] = 1;
        
        prev_node_id = node_id;
        node_id++;
      }
    }
#ifdef __PROFILE__
    mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif
#ifdef __DEBUG2__
    {
      fprintf(fout_db, "Adjacency lists:\n");
      Hash_UInt32_Int32::iterator it;
      int p;
      for(it=adjacency_list_pointer.begin(); it!=adjacency_list_pointer.end();
          it++){
        fprintf(fout_db, "label: %u\n", (*it).first);
        p = (*it).second;
        while(p!=-1){
          fprintf(fout_db, "neighbor: %u\n", 
                    adjacency_lists[p].neighbor_id);
          p = adjacency_lists[p].next_ptr;
        }
      }
    }
#endif //__DEBUG2__
  }
  
  // Make mergers. Save label mapping for f_threshold_seq
  plhs[0] = mxCreateCellMatrix(1, n_f_threshold);

#ifdef __DEBUG__
  fprintf(fout_db, "Start of mergers.\n");
  fflush(fout_db);
#endif //__DEBUG__

  {
    Q_Double_Label_Pair::iterator q_it, q_it1;
    double m, m1; // merge criterion for a pair
    unsigned int s0, s1; // segments to be merged
    unsigned int sl, sr, si;
    unsigned int n0, n1; // their respective neighbors
    int f_threshold_id=0;
    double f_threshold;
    unsigned int b, n; // boundary value and number of pixels
    int p, prev_p, next_p;
    Label_Pair lp_li, lp_il, lp_ir, lp_ri, lp_lr, lp_rl, lp_10, lp_01;
    q_it = merge_q.begin();
    m = (*q_it).first;
    tie(s0,s1) = id_2_label_pair((*q_it).second);
    is_valid_merge_q_pos[(*q_it).second]=0;
    merge_q.erase(q_it);
    lp_10 = label_pair_2_id(s1, s0);
    if(is_valid_merge_q_pos[lp_10]==1){
      merge_q.erase(merge_q_pos[lp_10]);
      is_valid_merge_q_pos[lp_10]=0;
    }
    while(f_threshold_id<n_f_threshold){
      f_threshold = f_thresholds[f_threshold_id];
      mexPrintf("f_threshold:%g, m:%g\n", f_threshold, m);
#ifdef __PROFILE__
      mexCallMATLAB(0, NULL, 0, NULL, "tic");
#endif
      while(q_it!=merge_q.end() && m<f_threshold){ // merge segments s0 and s1

        lp_01 = label_pair_2_id(s0, s1);
        if(label_adjacency[lp_01]==0 ||
           superpixel_sets.get_adam(s0)==superpixel_sets.get_adam(s1)){
          q_it = merge_q.begin();
          if(q_it==merge_q.end())
            break;
          m = (*q_it).first;
          tie(s0,s1) = id_2_label_pair((*q_it).second);
          is_valid_merge_q_pos[(*q_it).second]=0;
          merge_q.erase(q_it);
          lp_10 = label_pair_2_id(s1, s0);
          if(is_valid_merge_q_pos[lp_10]==1){
            merge_q.erase(merge_q_pos[lp_10]);
            is_valid_merge_q_pos[lp_10]=0;
          }
          continue; // already merged
        }
        
        // merge the superpixel sets of s1 into s0
        if(n_neighbor[s0]<n_neighbor[s1]){
          sl = s0;
          sr = s1;
        }
        else{
          sl = s1;
          sr = s0;
        }
        
        lp_lr = label_pair_2_id(sl, sr);
        lp_rl = label_pair_2_id(sr, sl);
        
        // mark this edge for deletion
        label_adjacency[lp_lr]=0; 
        label_adjacency[lp_rl]=0; 

#ifdef __DEBUG__
        mexPrintf("sl:%u, sr:%u, m:%g\n", sl, sr, m);
        fprintf(fout_db, "sl:%u, sr:%u, m:%g\n", sl, sr, m);
        fflush(fout_db);
#endif //__DEBUG__
        
        p = adjacency_list_pointer[sl];
        while(p>0){ // add the neighbors of sl to sr
          next_p = adjacency_lists[p].next_ptr;
          si = adjacency_lists[p].neighbor_id;
          lp_li = label_pair_2_id(sl,si);
          
#ifdef __DEBUG__
          fprintf(fout_db, "p:%d, next_p:%d, si:%u\n", p, next_p, si);
          fflush(fout_db);
#endif //__DEBUG__
          
          if(si==sr || label_adjacency[lp_li]==0 || n_boundary_pixel[lp_li]==0){
            // neighbor marked for deletion - delete and skip.
            p = next_p;
            continue;
          }
          
          lp_ir = label_pair_2_id(si, sr);
          if(label_adjacency[lp_ir]==0){ // si is not a neighbor of sr

#ifdef __DEBUG__
            fprintf(fout_db, "si is not a neighbor of sr ...\n");
            fflush(fout_db);
#endif //__DEBUG__
            
            lp_il = label_pair_2_id(si, sl);
            lp_ri = label_pair_2_id(sr, si);

            // add si to sr's neighbors
            label_adjacency[lp_ir] = 1;
            label_adjacency[lp_ri] = 1;
            
            n_boundary_pixel[lp_ir] = n_boundary_pixel[lp_il];
            n_boundary_pixel[lp_ri] = n_boundary_pixel[lp_li];
            sum_boundary_value[lp_ir] = sum_boundary_value[lp_il];
            sum_boundary_value[lp_ri] = sum_boundary_value[lp_li];

            if(is_valid_merge_q_pos[lp_il]==1){
              merge_q_pos[lp_ir] = merge_q_pos[lp_il];
              (*merge_q_pos[lp_ir]).second = lp_ir;
            }
            if(is_valid_merge_q_pos[lp_li]==1){            
              merge_q_pos[lp_ri] = merge_q_pos[lp_li];
              (*merge_q_pos[lp_ri]).second = lp_ri;
            }
            is_valid_merge_q_pos[lp_ir] = is_valid_merge_q_pos[lp_il];
            is_valid_merge_q_pos[lp_ri] = is_valid_merge_q_pos[lp_li];          
            is_valid_merge_q_pos[lp_il] = 0;
            is_valid_merge_q_pos[lp_li] = 0;

            label_adjacency[lp_li] = 0;
            label_adjacency[lp_il] = 0;

            adjacency_lists[p].next_ptr = adjacency_list_pointer[sr];
            adjacency_list_pointer[sr] = p;
            p = next_p;
            
#ifdef __DEBUG__
            fprintf(fout_db, "added si to sr.\n");
            fflush(fout_db);
#endif //__DEBUG__

          }
          else{ // si is a neighbor of sl and sr

#ifdef __DEBUG__
            fprintf(fout_db, "si is a neighbor of sr ... \n");
            fflush(fout_db);
#endif //__DEBUG__
            
            lp_il = label_pair_2_id(si, sl);
            lp_ri = label_pair_2_id(sr, si);
            
#ifdef __DEBUG__
            fprintf(fout_db, "combining mean boundary values ... ");
            fflush(fout_db);
#endif //__DEBUG__

            n_boundary_pixel[lp_ir] += n_boundary_pixel[lp_il];
            n_boundary_pixel[lp_ri] = n_boundary_pixel[lp_ir];
            
            sum_boundary_value[lp_ir] += sum_boundary_value[lp_il];
            sum_boundary_value[lp_ri] = sum_boundary_value[lp_ir];

#ifdef __DEBUG__
            fprintf(fout_db, "done.\n");
            fflush(fout_db);
#endif //__DEBUG__
           
            if(is_valid_merge_q_pos[lp_il]==1){
              merge_q.erase(merge_q_pos[lp_il]);
              is_valid_merge_q_pos[lp_il]=0;
            }
            if(is_valid_merge_q_pos[lp_li]==1){
              merge_q.erase(merge_q_pos[lp_li]);
              is_valid_merge_q_pos[lp_li]=0;
            }
            if(is_valid_merge_q_pos[lp_ir]==1){
              merge_q.erase(merge_q_pos[lp_ir]);
              is_valid_merge_q_pos[lp_ir]=0;
            }
            if(is_valid_merge_q_pos[lp_ri]==1){
              merge_q.erase(merge_q_pos[lp_ri]);
              is_valid_merge_q_pos[lp_ri]=0;
            }
                                                
#ifdef __DEBUG__
            fprintf(fout_db, "done.\n");
            fflush(fout_db);
#endif //__DEBUG__
            
            // insert in merge queue
            {

#ifdef __DEBUG__
              fprintf(fout_db, "inserting updated boundary into merge q ... ");
              fflush(fout_db);
#endif //__DEBUG__

              double m1 = (double)sum_boundary_value[lp_ir]/
                (double)n_boundary_pixel[lp_ir];
              Q_Double_Label_Pair::iterator q_it1;
              q_it1 = merge_q.find(m1); // keep priority keys unique.
              while(q_it1!=merge_q.end()){
                m1 += 0.00002*((double)(rand()%1000)); // epsilon increments
                q_it1 = merge_q.find(m1);
              }
              {
                const std::pair<double, Label_Pair> v(m1, lp_ir);
                bool is_inserted;
                tie(q_it1, is_inserted)  = merge_q.insert(v);
              }
              merge_q_pos[lp_ir] = q_it1;
              is_valid_merge_q_pos[lp_ir] = 1;

              m1 = (double)sum_boundary_value[lp_ri]/
                (double)n_boundary_pixel[lp_ri];
              // insert in merge queue
              q_it1 = merge_q.find(m1); // keep priority keys unique.
              while(q_it1!=merge_q.end()){
                m1 += 0.00002*((double)(rand()%1000)); // epsilon increments
                q_it1 = merge_q.find(m1);
              }
              {
                const std::pair<double, Label_Pair> v(m1, lp_ri);
                bool is_inserted;
                tie(q_it1, is_inserted)  = merge_q.insert(v);
              }
              merge_q_pos[lp_ri] = q_it1;
              is_valid_merge_q_pos[lp_ri] = 1;
              
#ifdef __DEBUG__
              fprintf(fout_db, "done.\n");
              fflush(fout_db);
#endif //__DEBUG__

            }
            p = next_p;

#ifdef __DEBUG__
            fprintf(fout_db, "merged info si-sl into sr.\n");
            fflush(fout_db);
#endif //__DEBUG__

          }
        }
        
        n_neighbor[sr] += n_neighbor[sl];
        superpixel_sets.merge(sl, sr);
       
        q_it = merge_q.begin();
        if(q_it==merge_q.end())
          break;
        m = (*q_it).first;
        tie(s0,s1) = id_2_label_pair((*q_it).second);
        is_valid_merge_q_pos[(*q_it).second]=0;
        merge_q.erase(q_it);
        lp_10 = label_pair_2_id(s1, s0);
        if(is_valid_merge_q_pos[lp_10]==1){
          merge_q.erase(merge_q_pos[lp_10]);
          is_valid_merge_q_pos[lp_10]=0;
        }
        
        f_threshold = f_thresholds[f_threshold_id];
      }
      
      // copy the current superpixel-to-segment mapping into
      // output variable
//      superpixel_sets.update_adams();
      mxArray * sp_2_seg_map_mx = mxCreateNumericMatrix(
        1, max_superpixel_id+1, mxUINT32_CLASS, mxREAL);
      {
        int i;
        unsigned int * sp_2_seg_map = (unsigned int*) mxGetPr(sp_2_seg_map_mx);
        for(i=0; i<max_superpixel_id+1; i++)
          sp_2_seg_map[i] = superpixel_sets.get_adam(i);
      }
      mxSetCell(plhs[0], f_threshold_id, sp_2_seg_map_mx);
      
#ifdef __PROFILE__
      mexCallMATLAB(0, NULL, 0, NULL, "toc");
#endif
      f_threshold_id++;
    }
  }

#ifdef __DEBUG__
  fprintf(fout_db, "STOP: segmentation_from_superpixels_mean_boundary_value_RAG2\n");
  fclose(fout_db);
#endif //__DEBUG__
 
  delete [] adjacency_lists;
  mexPrintf("STOP: segmentation_from_superpixels_mean_boundary_value_RAG2\n");
  return;
}
