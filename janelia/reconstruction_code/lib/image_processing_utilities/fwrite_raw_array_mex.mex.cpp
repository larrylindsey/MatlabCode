// Write a n_dimension array of elements of size_element bytes each to
// a file named file_name from a buffer pointed to by data.
// Wrapper for iput::io::fwrite_raw_array()
//
// Shiv N. Vitaladevuni
// Janelia Farm Research Campus, HHMI
// vitaladevunis@janelia.hhmi.org
//
// v0	06142009	init. code
//

#include <mex.h>

#include <file_input_output.h>

#define MIN(A,B)	((A)<(B)?(A):(B))
#define MAX(A,B)	((A)>(B)?(A):(B))

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){

  if(nrhs==0){
    mexPrintf("Usage err = fwrite_raw_array_mex(file_name, x)\n");
    mexPrintf("input params:\n");
    mexPrintf("\t1. file name (char array)\n");
    mexPrintf("\t2. variable to be written to file\n");
    mexPrintf("output:\n");
    mexPrintf("\t1. error code.\n");
    return;
  }
  if(nrhs!=2){
    mexErrMsgTxt("Wrong number of inputs\n");
    return;
  }
  if(nlhs>1){
    mexErrMsgTxt("Wrong number of outputs\n");
    return;
  }

  const mxArray * file_name_mx = prhs[0];
  const mxArray * variable_mx = prhs[1];
  
  unsigned int type_id;
  unsigned int size_element;
  { // get type of the variable
    int mx_type = mxGetClassID(variable_mx);
    if(mx_type == mxINT8_CLASS){
      type_id = iput::io::type_id_char();
      size_element = sizeof(char);
    }
    if(mx_type == mxUINT8_CLASS){
      type_id = iput::io::type_id_unsigned_char();
      size_element = sizeof(unsigned char);
    }
    if(mx_type == mxINT16_CLASS){
      type_id = iput::io::type_id_short_int();
      size_element = sizeof(short int);
    }
    if(mx_type == mxUINT16_CLASS){
      type_id = iput::io::type_id_unsigned_short_int();
      size_element = sizeof(unsigned short int);
    }
    if(mx_type == mxINT32_CLASS){
      type_id = iput::io::type_id_int();
      size_element = sizeof(int);
    }
    if(mx_type == mxUINT32_CLASS){
      type_id = iput::io::type_id_unsigned_int();
      size_element = sizeof(unsigned int);
    }
    if(mx_type == mxINT64_CLASS){
      type_id = iput::io::type_id_long_int();
      size_element = sizeof(long int);
    }
    if(mx_type == mxUINT64_CLASS){
      type_id = iput::io::type_id_unsigned_long_int();
      size_element = sizeof(unsigned long int);
    }
    if(mx_type == mxSINGLE_CLASS){
      type_id = iput::io::type_id_float();
      size_element = sizeof(float);
    }
    if(mx_type == mxDOUBLE_CLASS){
      type_id = iput::io::type_id_double();
      size_element = sizeof(double);
    }
  }
  iput::io::Type_Element type_element =
    iput::io::get_type_element(type_id, size_element);

  int n_dimension = mxGetNumberOfDimensions(variable_mx);
  int * dimensions = (int *) mxGetDimensions(variable_mx);

  char * file_name;
  {
    int l = mxGetN(file_name_mx)+1;
    file_name = (char *) new char[l];
    mxGetString(file_name_mx, file_name, l);
    mexPrintf("file_name: %s\n", file_name);
  }
  
  unsigned int err = iput::io::fwrite_raw_array(file_name, type_element,
                                                n_dimension, (int *) dimensions,
						(void*) mxGetPr(variable_mx));
  delete [] file_name;
  
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double * err_out = mxGetPr(plhs[0]);
  * err_out = (double) (err!=0);
  
  return;
}
