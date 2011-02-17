// Read a n_dimension array of elements of size_element bytes each
// from file named file_name. Memory is allocated and pointed to by
// data.  Wrapper for iput::io::fread_raw_array()
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	06142009	init. code
//

#include <string.h>
#include <mex.h>

#include <file_input_output.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

	if(nrhs==0){
		mexPrintf("Usage err = fread_raw_array_mex(file_name, x)\n");
		mexPrintf("input params:\n");
		mexPrintf("\t1. file name (char array)\n");
		mexPrintf("output:\n");
		mexPrintf("\t1. variable\n");
		mexPrintf("\t2. error code\n");
		return;
	}
	if(nrhs!=1){
		mexErrMsgTxt("Wrong number of inputs\n");
		return;
	}
	if(nlhs>2){
		mexErrMsgTxt("Wrong number of outputs\n");
		return;
	}

  const mxArray * file_name_mx = prhs[0];
  char * file_name;
  {
    int l = mxGetN(file_name_mx)+1;
    file_name = (char *) new char[l];
    mxGetString(file_name_mx, file_name, l);
  }
  
  iput::io::Type_Element type_element;
  int n_dimension;
  int * dimensions;
  void * buffer = iput::io::fread_raw_array(file_name, type_element,
                                            n_dimension, dimensions);

  delete [] file_name;
  
  bool error_flag = iput::io::is_error_type_element(type_element) ||
    buffer==NULL;
  if(nlhs>=2)
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);

  if(error_flag){
    plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    if(nlhs>=2)
      *(mxGetPr(plhs[1])) = error_flag;
    return;
  }

  unsigned int type_id = iput::io::get_type_id(type_element);
  unsigned int size_element = iput::io::get_size_element(type_element);
  mxClassID mx_type;
  { // get type of the variable
    if(type_id == iput::io::type_id_char())
      mx_type = mxINT8_CLASS;
    if(type_id == iput::io::type_id_unsigned_char())
      mx_type = mxUINT8_CLASS;
    if(type_id == iput::io::type_id_short_int())
      mx_type = mxINT16_CLASS;
    if(type_id == iput::io::type_id_unsigned_short_int())
      mx_type = mxUINT16_CLASS;
    if(type_id == iput::io::type_id_int())
      mx_type = mxINT32_CLASS;
    if(type_id == iput::io::type_id_unsigned_int())
      mx_type = mxUINT32_CLASS;
    if(type_id == iput::io::type_id_long_int())
      mx_type = mxINT64_CLASS;
    if(type_id == iput::io::type_id_unsigned_long_int())
      mx_type = mxUINT64_CLASS;
    if(type_id == iput::io::type_id_float())
      mx_type = mxSINGLE_CLASS;
    if(type_id == iput::io::type_id_double())
      mx_type = mxDOUBLE_CLASS;
  }
  
  plhs[0] = mxCreateNumericArray(n_dimension, (const int *) dimensions,
                                 mx_type, mxREAL);
  memcpy((void*) mxGetPr(plhs[0]), (void*) buffer,
         mxGetNumberOfElements(plhs[0])*size_element);

  if(nlhs>=2)
    *(mxGetPr(plhs[1])) = error_flag;

  {
    char * b = (char *) buffer;
    delete [] b;
  }
  return;
}
