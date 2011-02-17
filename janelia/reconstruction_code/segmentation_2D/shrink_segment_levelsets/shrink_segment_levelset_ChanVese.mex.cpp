// Shrink a given segmentation map to fit a boundary/interior probability map. This may be useful for modifying segmentation to improve linkage.
//
// Pratim Ghosh
// Summer Intern 2009,
// Chklovskii lab,
// Janelia Farm Research Campus,
// HHMI
//
// Shiv N. Vitaladevuni
// JFRC, HHMI
//
 

#include "stdio.h"
#include "math.h"
# include "stdlib.h"

#include "mex.h"  

#include <levelsets_ChanVese.h>
#include <hash_functions.h>
#include <queue>

# define eps 2.2204e-16
# define PI 3.1416

typedef std::pair<int, int> Coord_2D;

void get_distance_map_2D_4C(unsigned int * b_map, int * d_map, int width, int height){

  // initialize distance map
  {
    int x,y;
    for(x=0; x<width; x++){
      for(y=0; y<height; y++){
        d_map[y*width+x] = -1;
      }
    }
  }

  // initialize queue
  std::queue<Coord_2D> Q;
  {
    int x,y;
    for(x=0; x<width; x++){
      for(y=0; y<height; y++){
        if(b_map[y*width+x]==1){
          Q.push(Coord_2D(y,x));
          d_map[y*width+x] = 0;
        }
      }
    }
  }

  // start flood fill
  {
    while(!Q.empty()){
      Coord_2D p = Q.front();
      Q.pop();
      int y = p.first;
      int x = p.second;
      
      // east neighbor
      if(x>0){
        int y1 = y;
        int x1 = x-1;
        if(d_map[y1*width+x1]==-1){
          d_map[y1*width+x1] = d_map[y*width+x]+1;
          Q.push(Coord_2D(y1,x1));
        }
      }
      // west neighbor
      if(x<width-1){
        int y1 = y;
        int x1 = x+1;
        if(d_map[y1*width+x1]==-1){
          d_map[y1*width+x1] = d_map[y*width+x]+1;
          Q.push(Coord_2D(y1,x1));
        }
      }
      // north neighbor
      if(y>0){
        int y1 = y-1;
        int x1 = x;
        if(d_map[y1*width+x1]==-1){
          d_map[y1*width+x1] = d_map[y*width+x]+1;
          Q.push(Coord_2D(y1,x1));
        }
      }
      // south neighbor
      if(y<height-1){
        int y1 = y+1;
        int x1 = x;
        if(d_map[y1*width+x1]==-1){
          d_map[y1*width+x1] = d_map[y*width+x]+1;
          Q.push(Coord_2D(y1,x1));
        }
      }
    }
  }

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  
  if(nrhs==0){
    mexPrintf("Usage: segment_map_new = shrink_segment_levelset_ChanVese(...\n");
    mexPrintf("\tsegment_map_original, boundary_prob_map, lambda1, ...\n");
    mexPrintf("lambda2, nu, number_of_iter, c1, c2, delt);\n");
    mexPrintf("The default values for parameters:\n lambda1=lambda2=0.5, \n nu=0.5x255^2, \n number_of_iter=4000, \n c1=0.0, c2=255.0, \n delt=0.1. \n");
    return;
  }
  
  mexPrintf("START:  shrink_segment_levelset_ChanVese\n");

  double c1=0.0, c2=255.0, lambda1=0.5, lambda2=0.5;
  int  ITER=4000, h=1;
  double epsilon=1.0, nu=0.5*255*255, delt=0.1;
  double *labelMat, *Im;
  int r, c;

  labelMat = mxGetPr(prhs[0]);
  r = (int)mxGetM(prhs[0]);
  c = (int)mxGetN(prhs[0]);
  Im = mxGetPr(prhs[1]);
  if(nrhs>=3)
    lambda1 = mxGetScalar(prhs[2]);
  if(nrhs>=4)
    lambda2=mxGetScalar(prhs[3]);
  if(nrhs>=5)
    nu=mxGetScalar(prhs[4]);
  if(nrhs>=6)
    ITER=(int)mxGetScalar(prhs[5]);
  if(nrhs>=7)
    c1=(int)mxGetScalar(prhs[6]);
  if(nrhs>=8)
    c2=(int)mxGetScalar(prhs[7]);
  if(nrhs>=9)
    delt=mxGetScalar(prhs[8]);

  mexPrintf("lambda1: %g, lambda2: %g, nu: %g, ITER: %d, c1: %g, c2:%g\n",
            lambda1, lambda2, nu, ITER, c1, c2);
  mexPrintf("h: %d, epsilon: %g, delt: %g\n", h, epsilon, delt);
  mexEvalString("drawnow;");

  plhs[0] = mxCreateDoubleMatrix(r, c, mxREAL);
  double * Bigphi =  mxGetPr(plhs[0]);

  double *sublabelMat, *subIm, *returnlabelMat, *returnBiglabel, *returnIm, *init_phi, *phiNew, *phiTemp;
  int i, j, k, l, num, label, num_label, x, y, z, w, subr, subc, count, R, C;
  unsigned int *maskgrt0, *masklt0;
  int temp, *phigrt0, *philt0;
  
  // Get bounding boxes
  Hash_UInt32_UInt32 label_rmin, label_rmax, label_cmin, label_cmax;
  {
    for(i=0; i<r; i++) {
      for(j=0; j<c; j++) {
        num=(unsigned int)labelMat[(j*r)+i];
        if(label_rmin.find(num)==label_rmin.end())
          label_rmin[num]=r-1;
        
      
        if(label_cmin.find(num)==label_cmin.end())
          label_cmin[num]=c-1;

        label_rmin[num] = std::min(label_rmin[num],(unsigned int)i);
        label_rmax[num] = std::max(label_rmax[num],(unsigned int)i);
        label_cmin[num] = std::min(label_cmin[num],(unsigned int)j);
        label_cmax[num] = std::max(label_cmax[num],(unsigned int)j);
      }
    }
  }
  
  num_label=label_rmin.size();
  
  maskgrt0=new unsigned int [r*c];
  masklt0=new unsigned int [r*c];
  phigrt0=new int [r*c];
  philt0=new int [r*c];
  sublabelMat=new double [r*c];
  subIm=new double [r*c];
  init_phi=new double [r*c];
  phiNew=new double [r*c];
  phiTemp=new double [r*c];
  
  Hash_UInt32_UInt32::iterator hp;

  for(hp=label_rmin.begin();hp!=label_rmin.end();hp++) {
    label = (*hp).first;
    if (label!=0){
      //mexPrintf("label: %d\n", label);
      //mexEvalString("drawnow;");
      
      x=std::max((int)(label_rmin[label]-5), 0);
      y=std::min((int)(label_rmax[label]+5), r-1);
      
      z=std::max((int)(label_cmin[label]-5), 0);
      w=std::min((int)(label_cmax[label]+5), c-1);
      
      subr=y-x+1;
      subc=w-z+1;
      
      for (j=x; j<=y; j++) {
        for (k=z; k<=w; k++) {
          if ((int)labelMat[(k*c)+j]==label)
            sublabelMat[(k-z)*subr+(j-x)]=(double) label;
          else
            sublabelMat[(k-z)*subr+(j-x)]=0.0;
          subIm[(k-z)*subr+(j-x)]=Im[(k*c)+j];
        }
      }
      
      for (j=0; j<subr; j++){
        for (k=0; k<subc; k++){
          if (sublabelMat[k*subr+j]>0.0) {
            maskgrt0[k*subr+j]=1;
            masklt0[k*subr+j]=0;
          }
          else {
            maskgrt0[k*subr+j]=0;
            masklt0[k*subr+j]=1;
          }
        }
      }
      
      get_distance_map_2D_4C(maskgrt0, phigrt0, subr, subc);
      get_distance_map_2D_4C(masklt0, philt0, subr, subc);
      
      
      for (j=0; j<subr; j++){
        for (k=0; k<subc; k++){
          if (sublabelMat[k*subr+j]>0.0){
            init_phi[k*subr+j]=(double)philt0[k*subr+j]+0.5;
          }
          else {
            subIm[k*subr+j]=255.0;
            init_phi[k*subr+j]=-(double)phigrt0[k*subr+j]-0.5;
          }
        }
      }
      
      ChanVese(init_phi, phiNew, subIm, epsilon, nu, delt, subr,
               subc, c1, c2, lambda1, lambda2, h, ITER);

      for (j=x; j<=y; j++) {
        for (k=z; k<=w; k++) {
          if (phiNew[(k-z)*subr+(j-x)]>=0.0){
            if (label>Bigphi[(k*c)+j])
              Bigphi[(k*c)+j]=label;
          }
        }
      }
    }
  }
  
  delete [] sublabelMat;
  delete [] subIm;
  delete [] phiNew;
  delete [] phiTemp;
  delete [] init_phi;
  delete [] phigrt0;
  delete [] philt0;
  delete [] masklt0;
  delete [] maskgrt0;
  
  mexPrintf("STOP:  shrink_segment_levelset_ChanVese\n");
  
  return;
}
        
