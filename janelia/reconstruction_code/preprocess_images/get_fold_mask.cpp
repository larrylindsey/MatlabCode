/*
 * Find folds in EM images to prevent segmentation errors and for alignment
 *
 * Main routine by Lou Scheffer, Visiting Scientist, JFRC, HHMI, Nov. 25 2008.
 * MATLAB mex wrapper by Shiv Vitaladevuni, JFRC, HHMI.
 *
 * v0   11252008  init. code
 *
 */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <mex.h>

using namespace std;

typedef unsigned char uint8;

class MeanStd {  // for computing mean and standard deviation
  public:
    MeanStd(){sum = sum2 = 0.0; n=0;}
    void Element(double a){sum += a; sum2 += a*a; n++;}
    void Stats(double &avg, double &std){ avg = sum/n; std = sqrt(sum2/n - avg*avg);}
    int HowMany(){return n;}
    private:
      double sum, sum2;
      int n;
};


void FindFoldMask(uint8 *raster, int w, int h, uint8 *FoldMask)
{
  // Here, define the connected region(s) on the above layer.
  // First, find the mean and standard deviation of the 'real' (non-zero) pixels
  int npixels = w * h;
  MeanStd m;
  for(int i=0; i<npixels; i++) {
    if (raster[i] > 40)
      m.Element(raster[i]);
  }
  printf("%d real pixels, %f percent\n", m.HowMany(), m.HowMany()*100.0/npixels);
  double mean, std;
  m.Stats(mean, std);  // find existing statistics
  printf("Of the non-zero image points, mean= %f and std dev = %f\n", mean, std);
  
  // Make a normalized copy of the raster.  Do not copy any pixels next to black pixels
  // since they themselves are dubious.
  // Make 2 copies since the first (v) will be destroyed during connected region evaluation.
  vector<double> v;
  vector<double> vorig;
  for(int i=0; i<npixels; i++) {
    int y = i/w;
    int x = i-y*w;   // pixels next to background pixels should also be black
    int pix = raster[i];
    if (x-1 >= 0 && raster[i-1] == 0 ) pix = 0;
    if (x+1 <  w && raster[i+1] == 0 ) pix = 0;
    if (y-1 >= 0 && raster[i-w] == 0 ) pix = 0;
    if (y+1 <  h && raster[i+w] == 0 ) pix = 0;
    v.push_back((pix-mean)/std);
  }
  
  // Compute the threshold to use in determining folds.
  
  double thresh = 2.0;             // normally use -4 as a fold threshold
  if (mean - thresh*std < 1) {     // but if this is not possible, reduce value
    thresh = mean/std*0.95;
    printf("Forced to reduce threshold to %f\n", thresh);
  }
  
  // Now compute the fold mask
  
  for(int i=0; i<npixels; i++) {             // now a very simple thresholding
    FoldMask[i] = v[i] < -thresh ? 1 : 0;
  }
  
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
  
  const mxArray * image_mx = prhs[0];
	const int * size_0 = mxGetDimensions(image_mx);
	int width = size_0[1], height = size_0[0];
  unsigned char * image = (unsigned char *) mxGetPr(image_mx);
  
  mxArray * fold_mask_mx = mxCreateNumericMatrix(height, width, mxUINT8_CLASS, 0);
  unsigned char * fold_mask = (unsigned char *) mxGetPr(fold_mask_mx);
  
  FindFoldMask(image, height, width, fold_mask);
  
  plhs[0] = fold_mask_mx;
  return;
}

