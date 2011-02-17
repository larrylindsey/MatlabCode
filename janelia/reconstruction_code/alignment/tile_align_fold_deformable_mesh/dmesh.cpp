double DetailInitialThreshold = 0.10;      // When we call the detailed optimizer, it should be this good at least
double DetailFinalThreshold = 0.305;       // And at least this good when we finish
double TransAreClose = 24.0;           // two transforms are 'close' if within this value (at 2k scale)
double minl = 300;                     // we want edges between minl and 2*minl;
double MinTriangleArea = 10000.0;      // Minimum allowable size of a triangle (pixels)
double MinNormCorrelation = 0.34;      // Minimum normalized correlation to define a crude match
double MinNormGreat = 0.45;            // if we find one this good, we can stop searching
double MinMapArea = 50000;             // this much area should map from one region to another
double NumCorrPatchSize = 100000;      // Normalized corr not computed on sizes smaller than this
double ApproxScale = 1.0;              // relative size of images
double InitXScale = 1.0;               // separate scale for X, for processing tilts
double InitYScale = 1.0;               // separate scale for Y, for processing tilts
double InitialSkew = 0.0;	       // try this skew to start
double JustOne = 0.0;                  // generate just one transform
double OuterR = 0.0;		       // larger radius for difference of gaussians
double InnerR = 20.0;                  // inner radius for difference of gaussians
double InitialFourierMetric = 0.50;    // Fourier match after affine step
double FinalFourierMetric = 0.65;      // Final quality, averaged over all triangles
double AngleCenter = 0.0;              // Center of search range
double FinalYellowLimit = 0.500;        // Percentage of pixels that must be yellow in the final image

// these types defined in tiffio.h, but we don't want to create another
// dependency here.

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
#include "dmesh.h"
static FILE *of;

// The routine to set Parameters
void SetOne(const char *args, const char *name, double &rslt)
{
char str[256];
sprintf(str, "-%s=", name);
const char *p = strstr(args, str);
if (p == NULL)
    return;
p = strchr(p,'=');
rslt = atof(p+1);
}

void SetDMeshParams(const char *arg)
{
SetOne(arg, "DIT", DetailInitialThreshold);
SetOne(arg, "DFT", DetailFinalThreshold);
SetOne(arg, "TAC", TransAreClose);
SetOne(arg, "MNL", minl);
SetOne(arg, "MTA", MinTriangleArea);
SetOne(arg, "MNC", MinNormCorrelation);
SetOne(arg, "MNG", MinNormGreat);
SetOne(arg, "MMA", MinMapArea);
SetOne(arg, "NCP", NumCorrPatchSize);
SetOne(arg, "SCALE", ApproxScale);
SetOne(arg, "XSCALE", InitXScale);
SetOne(arg, "YSCALE", InitYScale);
SetOne(arg, "SKEW", InitialSkew);
SetOne(arg, "ONE", JustOne);
SetOne(arg, "OUTR", OuterR);
SetOne(arg, "INR", InnerR);
SetOne(arg, "IFM", InitialFourierMetric);
SetOne(arg, "FFM", FinalFourierMetric);
SetOne(arg, "CTR", AngleCenter);
SetOne(arg, "FYL", FinalYellowLimit);
}

//---------------------- End of correlation code -------------------------------------
void PrintTransform(FILE *of, tform &tr)
{
fprintf(of,"%7.4f %7.4f %8.2f   %7.4f %7.4f %8.2f\n", tr.t[0], tr.t[1], tr.t[2], tr.t[3], tr.t[4], tr.t[5]);
}

void InvertTrans(tform &inv, tform &b)
{
double det = b.t[0]*b.t[4] - b.t[1]*b.t[3];
inv.t[0] = b.t[4]/det; inv.t[1] = -b.t[1]/det; inv.t[2] = (b.t[5]*b.t[1]-b.t[4]*b.t[2])/det;
inv.t[3] =-b.t[3]/det; inv.t[4] =  b.t[0]/det; inv.t[5] = (b.t[2]*b.t[3]-b.t[0]*b.t[5])/det;
}

//Invert a 3x3 matrix.
double Invert3x3matrix( double i[3][3], double a[3][3])
{
double det = a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[1][0]*a[2][1]*a[0][2] -
             a[0][2]*a[1][1]*a[2][0] - a[0][1]*a[1][0]*a[2][2] - a[1][2]*a[2][1]*a[0][0];
i[0][0] =  (a[1][1]*a[2][2] - a[1][2]*a[2][1])/det;
i[1][0] = -(a[1][0]*a[2][2] - a[1][2]*a[2][0])/det;
i[2][0] =  (a[1][0]*a[2][1] - a[1][1]*a[2][0])/det;
i[0][1] = -(a[0][1]*a[2][2] - a[0][2]*a[2][1])/det;
i[1][1] =  (a[0][0]*a[2][2] - a[0][2]*a[2][0])/det;
i[2][1] = -(a[0][0]*a[2][1] - a[0][1]*a[2][0])/det;
i[0][2] =  (a[0][1]*a[1][2] - a[0][2]*a[1][1])/det;
i[1][2] = -(a[0][0]*a[1][2] - a[0][2]*a[1][0])/det;
i[2][2] =  (a[0][0]*a[1][1] - a[0][1]*a[1][0])/det;
return det;
}

// Compute r = a*b
void MultiplyTrans(tform &r, tform &a, tform &b)
{
r.t[0] = a.t[0]*b.t[0]+a.t[1]*b.t[3]; r.t[1] = a.t[0]*b.t[1]+a.t[1]*b.t[4]; r.t[2] = a.t[0]*b.t[2]+a.t[1]*b.t[5]+a.t[2];
r.t[3] = a.t[3]*b.t[0]+a.t[4]*b.t[3]; r.t[4] = a.t[3]*b.t[1]+a.t[4]*b.t[4]; r.t[5] = a.t[3]*b.t[2]+a.t[4]*b.t[5]+a.t[5];
}

void GetBaryValues(
 vector<Point> &cpts,         // control points
 vector<vector<double> > &lam, // point of interest, as multipliers of the control points.  If triangulated,
                              // at most 3 of these will be non-zero for each point
 vector<double> &raster,      // the raster we are mapping into
 int w,                       // the width of the raster
 vector<double> &spv,         // source pixel values.  Needed for computing derivatives
 vector<double> &v,           // the output pixel values
 vector<Point> &derivs)       // derivatives with respect to moving the control points
{
// initialize the derivatives.
int nctp = cpts.size();
derivs.resize(nctp);
for(int i=0; i<nctp; i++)
    derivs[i].x = derivs[i].y = 0.0;

int npts = lam.size();
v.resize(npts);
for(int i=0; i<npts; i++) {
    // find x, y by forming the linear combinations of the control points
    double x=0.0, y=0.0;
    for(int j=0; j<nctp; j++) {
        x += lam[i][j]*cpts[j].x;
        y += lam[i][j]*cpts[j].y;
	}
    // find the 4 points that surround the real valued coordinate
    if (x < 0.0 || x >= w-1 || y < 0.0 || y >= 4095) {  // bloot, should pass h
         v[i] = 0.0;  // anything outside the raster is 0.0
         continue;
         }
    int xll = int(x);
    int xur = xll+1;
    int yll = int(y);
    int yur = yll+1;
    // now do the interpolation
    double alpha = x-xll;
    double beta  = y-yll;
    double ll = raster[w*yll+xll];
    double ul = raster[w*yur+xll];
    double ur = raster[w*yur+xur];
    double lr = raster[w*yll+xur];
    double dx = lr-ll;
    double dy = ul-ll;
    double dxy = ll - lr - ul + ur;
    v[i] = ll + alpha*dx + beta*dy + alpha*beta*dxy;
    // now find the raw derivatives with respect to x and y
    double derx = dx +  beta*dxy;
    double dery = dy + alpha*dxy;
    // the derivatives of the correlation with respect to moving the control points
    for(int j=0; j<nctp; j++) {
        derivs[j].x += spv[i]*lam[i][j]*derx;
        derivs[j].y += spv[i]*lam[i][j]*dery;
	}
    }
}

void Transform(Point &p, tform &t)
{
double x = p.x*t.t[0] + p.y*t.t[1] + t.t[2];
double y = p.x*t.t[3] + p.y*t.t[4] + t.t[5];
p.x = x;
p.y = y;
}

void Transform(vector<Point> &v, tform &t)
{
for(int i=0; i<v.size(); i++) {
    Transform(v[i], t);
    }
}

// finds the correlation between two vectors, assuming each has mean 0 and std 1
// Since we are using this for correlation, just check the non-zero pixels
// Returns the number of non-zero pixels, if requested
double CorrNorm(vector<double> &a, vector<double> &b, int *nnz = NULL)
{
if (a.size() != b.size() ) {
    fprintf(of,"Sizes differ! %d %d\n", a.size(), b.size() );
    exit(42);
    }
double sum2 = 0.0;
int nz = 0;   // number of 0 entries
for(int i=0; i<a.size(); i++) {
    double prod = a[i]*b[i];
    sum2 += prod;
    if (fabs(prod)<1.0E-8)
	nz++;
    }
fprintf(of,"%d of %d were small\n", nz, a.size());
if (nnz != NULL)
    *nnz = a.size() - nz;
return sum2/(a.size() - nz);
}

// Finds statistics on a vector of doubles.
void Stats(vector<double> &v, double &avg, double &std)
{
MeanStd m;
int n = v.size();
for(int i=0; i<n; i++) 
    m.Element(v[i]);
m.Stats(avg, std);
}

// Normalizes a vector to have mean 0 and standard deviation 1
void Normalize(vector<double> &v)
{
double avg, std;
Stats(v, avg, std);
//fprintf(of,"Average %f, std dev %f\n", avg, std);
for(int i=0; i<v.size(); i++)
    v[i] = (v[i] - avg)/std;
}

//note the location of (x,y) in the array is different than the fftw doc, but this is OK in this
// case since we use only ffti(fft(d))
double LookFFT(double *fr, int N, int x, int y)
{
if (x < 0) x += N;
if (y < 0) y += N;
return fr[x + N*y]/N/N;
}

// same as above, but no divide by N
double LkFFT(double *fr, int N, int x, int y)
{
if (x < 0) x += N;
if (y < 0) y += N;
return fr[x + N*y];
}

// find where time partial image defined by pts and vals occurs in image2
// Look for the best one within radius RADIUS of (tx,ty)
double FindCorrelation(vector<Point> &pts, vector<double> &vals, vector<double> &image2,
 double &dx, double &dy, int tx, int ty, int radius)
{
// real to complex generates an array of size N *(n/2+1)
int N = 4096;
int M = N*(N/2+1);  // number of complex numbers in FFT of 2D real
fftw_plan p;
CD *out = (CD*) fftw_malloc(sizeof(CD) * M);
p = fftw_plan_dft_r2c_2d(N, N, (double *)(&image2[0]), (double (*)[2])out, FFTW_ESTIMATE);
fprintf(of,"start fft\n");
fftw_execute(p); // execute the plan
fprintf(of,"end fft\n");
fftw_destroy_plan(p);

// Now fft the partial image
vector<double> i1(N*N,0.0);
for(int i=0; i<pts.size(); i++) {
    int x = int(pts[i].x +0.5);  // Rounding should not be needed since integers of this size
    int y = int(pts[i].y +0.5);  // are exact in floating point
    i1[N*y+x] = vals[i];
    }
CD *out2 = (CD*) fftw_malloc(sizeof(CD) * M);
p = fftw_plan_dft_r2c_2d(N, N, (double *)(&i1[0]), (double (*)[2])out2, FFTW_ESTIMATE);
fprintf(of,"start fft #2\n");
fftw_execute(p); // execute the plan
fprintf(of,"end fft\n");
fftw_destroy_plan(p);

// Now multiply the original by the conjugate of the image.  Include filter
int yrow = N/2+1;
for(int i=0; i<M; i++) {
    int x = i/yrow;
    int y = i -yrow*x;
    if (x > N/2) x = x-N;
    double rad2 = (double(x)*double(x) + double(y)*double(y))/10000.0;  // radius squared
    double filt = 1.0;
    if (rad2 > 10)
        filt = 0.0;
    else
        filt = exp(-rad2);  // keep first 30 harmonics or so...
    out[i] = out[i] * conj(out2[i]) *filt;
    }

// reverse the FFT to find the lags
double *rslt = (double *) fftw_malloc(sizeof(double) * N*N);
p = fftw_plan_dft_c2r_2d(N, N, (double (*)[2])(&out[0]), rslt, FFTW_ESTIMATE);
fprintf(of,"start fft #3\n");
fftw_execute(p); // execute the plan
fprintf(of,"end fft\n");
fftw_destroy_plan(p);

// Now find the maximum value
double biggest = -1.0E30;
int bigx = -1;
int bigy = 0;
for(int i=0; i<N*N; i++) {
    int y = i/N;
    int x = i - y*N;
    if (x > N/2) x = x-N;
    if (y > N/2) y = y-N;
    if (abs(rslt[i])/N/N > 10.0*pts.size()) {  // too big to be true
         fprintf(of,"Very odd - %d %f %d\n", i, rslt[i]/N/N, pts.size());
         exit(43);
         }
    if (abs(x-tx) > radius || abs(y-ty) > radius)  // take only correlation within radius
	continue;
    if (rslt[i] > biggest) {
        biggest = rslt[i];
        bigx = x; bigy = y;
        }
    }
// for interest, write the landscape around the max
for(int iy=-6; iy <=6; iy += 2) {
    int ay = bigy + iy;
    fprintf(of,"biggest %f, y=%d, tx, ty, radius= %d %d %d\n", biggest, ay, tx, ty, radius);
    if (ay < 0)
        ay += N;
    for(int ix=-6; ix<=6; ix += 2) {
        int ax = bigx + ix;
        if (ax < 0)
            ax += N;
        fprintf(of,"%12.1f ",rslt[ax+ay*N]/N/N);
        }
    fprintf(of,"\n");
    }

fprintf(of,"Maximum correlation of %f at [%d,%d]\n", biggest/N/N/pts.size(), bigx, bigy);

// now interpolate to find a more accurate (we hope) peak.
double v = LookFFT(rslt, N, bigx, bigy);
double mi = LookFFT(rslt, N, bigx-1, bigy) - v;  // relative value at minus 1 coordinate
double pl = LookFFT(rslt, N, bigx+1, bigy) - v;  // same for +1 coordinate
double a = (mi+pl)/2;
double b = (pl-mi)/2;  // local equation now a*x^2 + b*x + c (but c == 0) by construction
double deltax = -b/(2*a);
fprintf(of,"%f %f %f %f %f -> dx = %f\n", mi, v, pl, a, b, deltax);
mi = LookFFT(rslt, N, bigx, bigy-1) - v;  // relative value at minus 1 coordinate
pl = LookFFT(rslt, N, bigx, bigy+1) - v;  // same for +1 coordinate
a = (mi+pl)/2;
b = (pl-mi)/2;  // local equation now a*x^2 + b*x + c (but c == 0) by construction
double deltay = -b/(2*a);
// Free the arrays we allocated
fftw_free(out);
fftw_free(out2);
fftw_free(rslt);
dx = bigx + deltax; dy = bigy + deltay;
return biggest/N/N/pts.size();
}

//Compute Normalized cross correlation straight from the definition.  Not efficient - only for debugging
double NormCrossDirect(vector<double> &a, vector<double> &b, int S,
  int ax1, int ay1, int ax2, int ay2, int bx1, int by1, int bx2, int by2)
{
int nx = ax2-ax1+1;
int ny = ay2-ay1+1;
double suma = 0.0, sumb = 0.0;
for(int x=0; x<nx; x++) {
    for(int y=0; y<ny; y++) {
        suma += a[ax1+x + (ay1+y)*S];
        sumb += b[bx1+x + (by1+y)*S];
        }
    }
double avga = suma/nx/ny;
double avgb = sumb/nx/ny;
double sumn=0.0, sumd1=0.0, sumd2=0.0;
for(int x=0; x<nx; x++) {
    for(int y=0; y<ny; y++) {
        double va = a[ax1+x + (ay1+y)*S];
        double vb = b[bx1+x + (by1+y)*S];
        sumn += (va-avga)*(vb-avgb);
        sumd1 += (va-avga)*(va-avga);
        sumd2 += (vb-avgb)*(vb-avgb);
        }
    }
double prod = sumd1*sumd2;
double rslt = prod < 1.0e-9 ? 0.0 : sumn/sqrt(prod);
fprintf(of,"\nNCD: s %f %f, a %f %f, ss %f %f %f, r %f\n", suma, sumb, avga, avgb, sumn, sumd1, sumd2, rslt);
return rslt;
}

//Given two images, each starting at 0,0.  Image 1 is WxH and image 2 is size W2xH2.  
// Assuming image 1 is shifted by adding (x,y) to its coordinates, then
// where does the overlap region fall in each image?
void BoxesFromShifts(int w, int h, int w2, int h2, int x, int y, 
 int &xl1, int &yb1, int &xr1, int &yt1, int &xl2, int &yb2, int &xr2, int &yt2)
{
// after shifting, bbox for image1 will be [x, x+w-1] in x, and [y,y+h-1] in y
xl2 = max(0, x);       // greater of left edges
xr2 = min(w2-1, x+w-1); // least of the rights
yb2 = max(0, y);       // greater of the bottoms
yt2 = min(h2-1, y+h-1); // lesser of the tops

xl1 = xl2-x;  // bounding box in image 1 coordinates is offset
xr1 = xr2-x;
yb1 = yb2-y;
yt1 = yt2-y;
}

// Look up in a cumulative table of doubles
double IntegralTable(vector<double> &t, int S, int x1, int y1, int x2, int y2)
{
double rslt = t[x2 + y2*S];
if (x1 > 0)
    rslt -= t[x1-1 + y2*S];
if (y1 > 0)
    rslt -= t[x2 + (y1-1)*S];
if (x1 > 0 && y1 > 0)
    rslt += t[(x1-1) + (y1-1)*S];
return rslt;
}

// same for ints
int IntegralTable(vector<int> &t, int S, int x1, int y1, int x2, int y2)
{
int rslt = t[x2 + y2*S];
if (x1 > 0)
    rslt -= t[x1-1 + y2*S];
if (y1 > 0)
    rslt -= t[x2 + (y1-1)*S];
if (x1 > 0 && y1 > 0)
    rslt += t[(x1-1) + (y1-1)*S];
return rslt;
}

// Compute the normalized cross-correlation.
// returns dx and dy.  These must be added to the coordinates of a point in image 1 to get them to match image2.
// Note: Last argument is a function that tells which matches are legal, in terms of x and y overlap size.
// pts and vals are describe the image to be correlated.  The other image is in a NxN vector of doubles.
// dx and dy are the returned offsets.
// tx, ty, and define a region.  We find the best match inside this region
// S is the size of the original image
// flog is the log file
// LegalRegion is a function that tells if this overlap (size sx by sy) should be considered.  'arg' is
//  also passed to this function.
// 'ftc' is a cache of the fourier transform of the second image.  If it is the right size, we assume
//   it has the right data.  To enforce re-computation, make it zero size
double FindNormCorrelation(vector<Point> &pts, vector<double> &vals, vector<Point> &ip2, vector<double> &iv2,
double &dx, double &dy, int tx, int ty, int radius, FILE *flog, 
bool (*LegalRegion)(int, int, void *), void *arg,   // function for checking if region dimensions are legal
bool (*LegalCounts)(int, int, void *), void *arg2,  // function for checking if point counts are legal
vector<CD> &ftc)
{
of = flog;

// Find the bounding box of the point list
int xmin = BIG, ymin = BIG, xmax = -BIG, ymax = -BIG;
for(int j=0; j<pts.size(); j++) {
    xmin = min(xmin, int(floor(pts[j].x)));
    xmax = max(xmax, int(ceil(pts[j].x)));
    ymin = min(ymin, int(floor(pts[j].y)));
    ymax = max(ymax, int(ceil(pts[j].y)));
    }
fprintf(of,"NormCorr: region size is [%d %d] in x, [%d %d] in y\n", xmin, xmax, ymin, ymax);
int w = xmax-xmin+1;
int h = ymax-ymin+1;  // size of region to be matched

// Find the bounding box of the image 2, the target
int xmin2 = BIG, ymin2 = BIG, xmax2 = -BIG, ymax2 = -BIG;
for(int j=0; j<ip2.size(); j++) {
    xmin2 = min(xmin2, int(floor(ip2[j].x)));
    xmax2 = max(xmax2, int(ceil (ip2[j].x)));
    ymin2 = min(ymin2, int(floor(ip2[j].y)));
    ymax2 = max(ymax2, int(ceil (ip2[j].y)));
    }
fprintf(of,"NormCorr: target image is [%d %d] in x, [%d %d] in y\n", xmin2, xmax2, ymin2, ymax2);
int w2 = xmax2-xmin2+1;
int h2 = ymax2-ymin2+1;  // size of region to be matched

// Find N, the size we need for the FFTs.  We will always use square FFTs
int N = 1;
for(; N < w+w2-1 || N < h+h2-1; )
    N = N*2;
fprintf(of, "N = %d\n", N);
int M = N*(N/2+1);  // number of complex numbers in FFT of 2D real

vector<double>i2(N*N,0.0);
// create the image2 from the point list.  Translate down to 0.0
// This is needed to create the running sums.  These could be cached as well.
for(int i=0; i<ip2.size(); i++) {
    int x = int(floor(ip2[i].x));
    double alpha = ip2[i].x - x;
    x -= xmin2;  
    int y = int(floor(ip2[i].y));
    double beta = ip2[i].y - y;
    y -= ymin2;  
    int f = N*y+x;   // first point, the lower left
    i2[f]     += (1-alpha)*(1-beta)*iv2[i];
    i2[f+1]   += (  alpha)*(1-beta)*iv2[i];  // next in X
    i2[f+N]   += (1-alpha)*(  beta)*iv2[i];  // next in Y
    i2[f+N+1] += (  alpha)*(  beta)*iv2[i];  // both up and over
    }

fftw_plan p;
if (ftc.size() != M) { // we have no cached FFT of image 2.  Compute one
    // real to complex generates an array of size N *(n/2+1)
    ftc.resize(M);
    //CD *out = (CD*) fftw_malloc(sizeof(CD) * M);
    p = fftw_plan_dft_r2c_2d(N, N, (double *)(&i2[0]), (double (*)[2])&(ftc[0]), FFTW_ESTIMATE);
    fprintf(of,"start fft\n");
    fftw_execute(p); // execute the plan
    fprintf(of,"end fft\n");
    fftw_destroy_plan(p);
    }

// create the image1 from the point list.  Translate down to 0.0
vector<double> i1(N*N,0.0);
for(int i=0; i<pts.size(); i++) {
    int x = int(floor(pts[i].x));
    double alpha = pts[i].x - x;
    x -= xmin;  
    int y = int(floor(pts[i].y));
    double beta = pts[i].y - y;
    y -= ymin;  
    int f = N*y+x;   // first point, the lower left
    i1[f]     += (1-alpha)*(1-beta)*vals[i];
    i1[f+1]   += (  alpha)*(1-beta)*vals[i];  // next in X
    i1[f+N]   += (1-alpha)*(  beta)*vals[i];  // next in Y
    i1[f+N+1] += (  alpha)*(  beta)*vals[i];  // both up and over
    }
// we need 4 arrays, each summarizing a value from (0,0) to (x,y)
// create 2 more, so we can tell the number of non-zero pixels in each portion
vector<double> i1sum (w*h,0.0);
vector<double> i1sum2(w*h,0.0);
vector<int>    i1nz(w*h, 0);
vector<double> i2sum (w2*h2,0.0);
vector<double> i2sum2(w2*h2,0.0);
vector<int>    i2nz(w2*h2,0);
i1sum[0]  = i1[0];        // first element
i1sum2[0] = i1[0]*i1[0];
i1nz[0] = (i1[0] != 0.0);
for(int i=1; i<w; i++) {
   i1sum[i]  = i1sum[i-1]  + i1[i];        // here x = i
   i1sum2[i] = i1sum2[i-1] + i1[i]*i1[i];
   i1nz[i] = i1nz[i-1] + (i1[i] != 0.0);
   }
for(int i=1; i<h; i++) {
   int j = i*w;  // stride of the array is w
   int k = i*N;  // but stride of the image is N
   i1sum[j]  = i1sum[j-w]  + i1[k];        // here y = i
   i1sum2[j] = i1sum2[j-w] + i1[k]*i1[k];
   i1nz[j] = i1nz[j-w] + (i1[k] != 0.0);
   }
i2sum[0]  = i2[0];
i2sum2[0] = i2[0]*i2[0];
i2nz[0] = (i2[0] != 0.0);
for(int i=1; i<w2; i++) {
   i2sum[i]  = i2sum[i-1]  + i2[i];        // here x = i
   i2sum2[i] = i2sum2[i-1] + i2[i]*i2[i];
   i2nz[i] = i2nz[i-1] + (i2[i] != 0.0);
   }
for(int i=1; i<h2; i++) {
   int j = i*w2;  // stride of the array is w
   int k = i*N;  // but stride of the image is N
   i2sum[j]  = i2sum[j-w2]  + i2[k];        // here y = i
   i2sum2[j] = i2sum2[j-w2] + i2[k]*i2[k];
   i2nz[j] = i2nz[j-w2] + (i2[k] != 0.0);
   }
for(int y=1; y<h; y++) {
    int i = y*w;  // address of first entry in the row
    int k = y*N;  // same for the image
    for(int x=1; x<w; x++) {
        int j = i+x; // index of value we want to set
        int m = k+x; // index in source array
	i1sum[j]  = i1sum[j-1]  + i1sum[j-w]  - i1sum[j-w-1]  + i1[m];
	i1sum2[j] = i1sum2[j-1] + i1sum2[j-w] - i1sum2[j-w-1] + i1[m]*i1[m];
        i1nz[j]   = i1nz[j-1]   + i1nz[j-w]   - i1nz[j-w-1] + (i1[m] != 0.0);
        }
    }
for(int y=1; y<h2; y++) {
    int i = y*w2;
    int k = y*N;
    for(int x=1; x<w2; x++) {
        int j = i+x; // index of value we want to set
        int m = k+x; // index in source array
	i2sum[j]  = i2sum[j-1]  + i2sum[j-w2]  - i2sum[j-w2-1]  + i2[m];
	i2sum2[j] = i2sum2[j-1] + i2sum2[j-w2] - i2sum2[j-w2-1] + i2[m]*i2[m];
        i2nz[j]   = i2nz[j-1]   + i2nz[j-w2]   - i2nz[j-w2-1]   + (i2[m] != 0.0);
        }
    }
//for(int y=h-1; y>= 0; y--) {
    //for(int x=0; x < w; x++) {
	//fprintf(of,"%7.0f", i1[y*N+x]);
        //}
    //fprintf(of,"\n");
    //}
//for(int y=h-1; y>= 0; y--) {
    //for(int x=0; x < w; x++) {
	//fprintf(of,"%7.0f", i1sum2[y*w+x]);
        //}
    //fprintf(of,"\n");
    //}
//for(int y=min(8,h2-1); y>= 0; y--) {
    //for(int x=0; x <= 8 && x < w2; x++) {
	//fprintf(of,"%7.0f", i2sum2[y*w2+x]);
        //}
    //fprintf(of,"\n");
    //}
// Now fft the partial image
CD *out2 = (CD*) fftw_malloc(sizeof(CD) * M);
p = fftw_plan_dft_r2c_2d(N, N, (double *)(&i1[0]), (double (*)[2])out2, FFTW_ESTIMATE);
fprintf(of,"start fft #2\n");
fftw_execute(p); // execute the plan
fprintf(of,"end fft\n");
fftw_destroy_plan(p);

// Now multiply the original by the conjugate of the image.
int yrow = N/2+1;
for(int i=0; i<M; i++) {
    int x = i/yrow;
    int y = i -yrow*x;
    if (x > N/2) x = x-N;
    out2[i] = ftc[i] * conj(out2[i]);
    }

// reverse the FFT to find the lags
double *rslt = (double *) fftw_malloc(sizeof(double) * N*N);
p = fftw_plan_dft_c2r_2d(N, N, (double (*)[2])(&out2[0]), rslt, FFTW_ESTIMATE);
fprintf(of,"start fft #3\n");
fftw_execute(p); // execute the plan
fprintf(of,"end fft\n");
fftw_destroy_plan(p);

// Now find the maximum value.  We'll also keep an array by size;
// We index into the array by int(log(N)*10), where N is the number of overlap pixels.
int lmax = (int)(10.0*log(double(N*N)));  // a very conservative upper bound
vector<double>max_by_size(lmax+1, 0.0);

// Now for the conventional biggest
double biggest = -1.0E30;
int bigx = -1;
int bigy = 0;
for(int i=0; i<N*N; i++) {
    int y = i/N;
    int x = i - y*N;
    if (x > N/2) x = x-N;
    if (y > N/2) y = y-N;
    if (x <= -w || x >= w2 || y <= -h || y >= h2) {
        rslt[i] = 0.0;
        continue;
	}
    // for each shift, find the overlap region in the coordinates of each image
    int xl1, yb1, xl2, yb2;
    int xr1, yt1, xr2, yt2;
    BoxesFromShifts(w, h, w2, h2, x, y, xl1, yb1, xr1, yt1, xl2, yb2, xr2, yt2);

    // Compute the cross-product from the definition, for debugging only
    //NormCrossDirect(i1, i2, 4096, xl1, yb1, xr1, yt1, xl2, yb2, xr2, yt2);
 
    //See if it's a legal overlap region.
    int sx = xr1-xl1+1;  // size in x
    int sy = yt1-yb1+1;  // size in y
    if (!(*LegalRegion)(sx, sy, arg)) {
        rslt[i] = 0.0;
	continue;
	}

    // See if there are enough of each kind of points in each region
    int i1count    = IntegralTable(i1nz,   w, xl1, yb1, xr1, yt1);
    int i2count    = IntegralTable(i2nz,   w2, xl2, yb2, xr2, yt2);
    if (LegalCounts != NULL && !(*LegalCounts)(i1count, i2count, arg2)) {
	rslt[i] = 0.0;
	continue;
	}
    // OK, region is legal.  Find the normalized cross correlation
    double n = sx*sy;
    double im1sum  = IntegralTable(i1sum,  w, xl1, yb1, xr1, yt1);
    double im1sum2 = IntegralTable(i1sum2, w, xl1, yb1, xr1, yt1);
    double im2sum  = IntegralTable(i2sum,  w2, xl2, yb2, xr2, yt2);
    double im2sum2 = IntegralTable(i2sum2, w2, xl2, yb2, xr2, yt2);
    double num = im1sum*im2sum/n -im1sum/n*im2sum -im2sum/n*im1sum + rslt[i]/N/N;
    double d1 = im1sum2 - im1sum/n*im1sum;
    double d2 = im2sum2 - im2sum/n*im2sum;
    double d = d1*d2;
    double r = d < 1.0E-9 ? 0.0: num/sqrt(d);  // if the denominator is very small, probable have underlying constant
    rslt[i] = r;
    //printf("final result = %f\n", r);
    //fprintf(of,"shift %d %d, i1 (%d %d) to (%d %d), i2 (%d %d) to (%d %d)\n",
     //x,y, xl1,yb1, xr1,yt1,  xl2,yb2, xr2,yt2);
        
    if (abs(r) > 1.0001) {  // too big to be true
         fprintf(of,"Very odd - i=%d rslt[i]=%f pts.size()=%d\n", i, rslt[i], pts.size());
	fprintf(of,"shift %d %d, i1 (%d %d) to (%d %d), i2 (%d %d) to (%d %d)\n",
	 x,y, xl1,yb1, xr1,yt1,  xl2,yb2, xr2,yt2);
        fprintf(of, "sums %f %f %f %f\nnum %f, d1 d2 %f %f\n", im1sum, im1sum2, im2sum, im2sum2, num, d1, d2);
	NormCrossDirect(i1, i2, 4096, xl1, yb1, xr1, yt1, xl2, yb2, xr2, yt2);
         exit(44);
         }
    if (abs(x-tx) > radius || abs(y-ty) > radius)  // take only correlation within radius
	continue;
    if (rslt[i] > biggest) {
        biggest = rslt[i];
        bigx = x; bigy = y;
        }
    int np = max(1, min(i1count, i2count));  // take the min of the point counts, but always >= 1
    int indx = (int)(10.0*log(double(np)));
    max_by_size[indx] = fmax(max_by_size[indx], rslt[i]);
    }
if (biggest < -2.0) {
    fprintf(of, "No legal subregions at all...\n");
    fprintf(of,"Maximum correlation of %f at [%d,%d]\n", 0.0, 0, 0);  // print so greps are helpful
    fftw_free(out2);
    fftw_free(rslt);
    return 0.0;
    }
// for interest, write the landscape around the max
for(int iy=-3; iy <=3; iy += 1) {
    int ay = bigy + iy;
    //fprintf(of,"biggest %f, y=%d, tx, ty, radius= %d %d %d\n", biggest, ay, tx, ty, radius);
    if (ay < 0)
        ay += N;
    for(int ix=-3; ix<=3; ix += 1) {
        int ax = bigx + ix;
        if (ax < 0)
            ax += N;
        fprintf(of,"%12.6f ",rslt[ax+ay*N]);
        }
    fprintf(of,"\n");
    }
// For debugging, print the results of correlation by region size.
// Note that a smaller region with less correlation is not a potential match,
// but if we find a larger region with a correlation almost as good, that's potentially
// a better match.
double bc = -10.0; // biggest correlation
int we = 0;        // which entry had it?
for(int i=0; i<max_by_size.size(); i++) {
    if (max_by_size[i] > bc) {
	bc = max_by_size[i];
        we = i;
	}
    }
// Now print the entries with bigger size and comparable correlation
// Mpy by 64 to get into pixel counts in the 2K working image.
const double PT = 0.8;  // Print threshold
for(int i=we+1; i<max_by_size.size(); i++) {
    if (max_by_size[i]/bc >= PT) {
        int i1 = (int) ceil(64.0*exp(i/10.0));
        int i2 = (int)floor(64.0*exp((i+1)/10.0));
	fprintf(of, "Possible bigger area match: %8d - %8d : %8.2f\n", i1, i2, max_by_size[i]);
	}
    }

fprintf(of,"Maximum correlation of %f at [%d,%d]\n", biggest, bigx, bigy);
//double limit = 0.9*biggest;
//for(int x=-w+1; x<w2-1; x++) {
    //for(int y=-h+1; y<h2-1; y++) {
	//double a = LkFFT(rslt, N, x, y);
	//if (a > limit && a > LkFFT(rslt, N, x-1, y) && a > LkFFT(rslt, N, x+1, y) &&
            //a > LkFFT(rslt, N, x, y-1) && a > LkFFT(rslt, N, x, y+1) ){
	    //fprintf(of, "Local max at %d %d, value %f\n", x, y, a);
	    //}
	//}
    //}

// now interpolate to find a more accurate (we hope) peak.
double v = LookFFT(rslt, N, bigx, bigy);
double mi = LookFFT(rslt, N, bigx-1, bigy) - v;  // relative value at minus 1 coordinate
double pl = LookFFT(rslt, N, bigx+1, bigy) - v;  // same for +1 coordinate
double a = (mi+pl)/2;
double b = (pl-mi)/2;  // local equation now a*x^2 + b*x + c (but c == 0) by construction
double deltax = -b/(2*a);
fprintf(of,"%f %f %f %f %f -> dx = %f\n", mi, v, pl, a, b, deltax);
mi = LookFFT(rslt, N, bigx, bigy-1) - v;  // relative value at minus 1 coordinate
pl = LookFFT(rslt, N, bigx, bigy+1) - v;  // same for +1 coordinate
a = (mi+pl)/2;
b = (pl-mi)/2;  // local equation now a*x^2 + b*x + c (but c == 0) by construction
double deltay = -b/(2*a);
// Free the arrays we allocated
fftw_free(out2);
fftw_free(rslt);
dx = bigx + deltax - xmin + xmin2; dy = bigy + deltay - ymin + ymin2;
return biggest;
}

int num = 1;  // write out converted images as 1.tif, 2.tif, etc.


// Improve the correlation by tweaking the location of the control points.
// returns the best correlation obtained.
double ImproveControlPts(
 vector<Point> &cpts,              // the control points
 vector<vector<double> > &lambda,  // pixels corrds expressed as linear combo of control points
 vector<double> &spv,              // the values at these pixels
 vector<double> &image2,           // the image we are mapping to
 int w,                            // and its width
 FILE *fsum,                       // Most important stuff goes here
 const char *describe,             // string to describe operation
 double threshold)                 // just for printing, not used otherwise
{
fprintf(of,"Contains %d pixels\n", spv.size() );
Normalize(spv);
   
double best_so_far = 0.0;

vector<double> Tpixels;          // Computed target pixel values go here
vector<Point> cderivs;           // and control point derivatives go here
GetBaryValues(cpts, lambda, image2, 4096, spv, Tpixels, cderivs);
int nnz;
double corr = CorrNorm(spv, Tpixels, &nnz);
//now do an initial check for plausibility
if (corr < DetailInitialThreshold) {  // if unexpected, dump out data passed in for debugging
    fprintf(of,"STAT:Correlation less than %f at the start of iterations, was %f\n", DetailInitialThreshold, corr);
    fprintf(fsum,"STAT:Correlation less than %f at the start of iterations, was %f\n", DetailInitialThreshold, corr);
    fprintf(of,"Control points are:");
    for(int i=0; i<cpts.size(); i++)
        fprintf(of,"(%.3f %.3f) ", cpts[i].x, cpts[i].y);
    fprintf(of,"corr=%f\n", corr);
    //for(int i=0; i<lambda.size(); i++) {
	//fprintf(of,"---i=%d\n",i);
	//for(int j=0; j<lambda[i].size(); j++)
	    //fprintf(of,"%.4f ", lambda[i][j]);
        //fprintf(of,"\n spv=%f Tpixels=%f\n", spv[i], Tpixels[i]);
        //}
    return 0.0;
    }
fprintf(of,"STAT: Initial %s correlation %f\n", describe, corr);
fprintf(fsum,"STAT: Initial %s correlation %f\n", describe, corr);

// Now try to tweak the control points for a good match
for(double step=10; step > 0.05; ) {
    fprintf(of,"\n");
    fprintf(of,"Control points are:");
    for(int i=0; i<cpts.size(); i++)
        fprintf(of,"(%f %f) ", cpts[i].x, cpts[i].y);
    fprintf(of,"corr=%f, step is %f\n", corr, step);
    // compute a downhill vector of length 'step'.
    double sum2 = 0.0;  // sum of squares
    for(int i=0; i<cderivs.size(); i++) {
        Point p = cderivs[i];
	sum2 += p.x*p.x + p.y*p.y;
        }
    double mag = sqrt(sum2);  // length of the vector
    // find the new (and hopefully improved) control points
    vector<Point> newpts = cpts;
    for(int i=0; i<cpts.size(); i++) {
       newpts[i].x += step*cderivs[i].x/mag;
       newpts[i].y += step*cderivs[i].y/mag;
       }
    // is the new spot better?
    vector<Point>  newderivs;
    GetBaryValues(newpts, lambda, image2, 4096, spv, Tpixels, newderivs);
    double c = CorrNorm(spv, Tpixels, &nnz);
    if (c > corr) { // new point is better
	corr = c;
        cpts = newpts;
        cderivs = newderivs;
        }
    else
	step = step / 2;
    }
fprintf(of,   "STAT: Final %s correlation %f, (threshold %f)\n", describe, corr, threshold);
fprintf(fsum, "STAT: Final %s correlation %f, (threshold %f)\n", describe, corr, threshold);
return corr;
}

void PrintMatrix(double a[3][3])
{
for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++)
	fprintf(of,"%12.6f ", a[i][j]);
    fprintf(of,"\n");
    }
}

void TryNewOptimizer(vector<Point> &plist, vector<double> &spv, vector<double> &image2, tform &t, FILE *flog)
{
// compute the bounding box
fprintf(of,"\n---------- Try new optimizer on %d points----------------\n", plist.size());
double xmin = BIG, ymin = BIG, xmax = -BIG, ymax = -BIG;
for(int j=0; j<plist.size(); j++) {
    xmin = min(xmin, plist[j].x);
    xmax = max(xmax, plist[j].x);
    ymin = min(ymin, plist[j].y);
    ymax = max(ymax, plist[j].y);
    }
fprintf(of,"region size is [%f %f] in x, [%f %f] in y\n", xmin, xmax, ymin, ymax);

// create 3 control points
vector<Point> cpts;
cpts.push_back(Point(xmin, ymin));
cpts.push_back(Point(xmax, ymin));
cpts.push_back(Point((xmin+xmax)/2, ymax));
fprintf(of,"Control points are (%f %f) (%f %f) (%f %f)\n", cpts[0].x, cpts[0].y, cpts[1].x, cpts[1].y, cpts[2].x, cpts[2].y);

// find each point as a linear combination of control points
double a[3][3];
a[0][0] = cpts[0].x; a[0][1] = cpts[1].x; a[0][2] = cpts[2].x;
a[1][0] = cpts[0].y; a[1][1] = cpts[1].y; a[1][2] = cpts[2].y;
a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
double inv[3][3];
PrintMatrix(a);
Invert3x3matrix(inv, a);
PrintMatrix(inv);
vector<vector<double> > lambda;
for(int j=0; j<plist.size(); j++) {
    //fprintf(of," Point is (%f %f)\n", plist[j].x, plist[j].y);
    vector<double> lam;
    lam.push_back(inv[0][0]*plist[j].x + inv[0][1]*plist[j].y + inv[0][2]*1.0);
    lam.push_back(inv[1][0]*plist[j].x + inv[1][1]*plist[j].y + inv[1][2]*1.0);
    lam.push_back(inv[2][0]*plist[j].x + inv[2][1]*plist[j].y + inv[2][2]*1.0);
    //fprintf(of," lambdas are %f %f %f\n", lam[0], lam[1], lam[2]);
    lambda.push_back(lam);
    }

// Transform the control points to the target frame
vector<Point> orig = cpts;
Transform(cpts, t);

// call the optimizer
ImproveControlPts(cpts, lambda, spv, image2, 4096, flog, "new opt", 1.0);

// Now, find a transformation that maps ORIG into the new cpts
// first, create a transform that maps a unit right triangle to the original pts
tform o(orig[1].x-orig[0].x, orig[2].x-orig[0].x, orig[0].x,
        orig[1].y-orig[0].y, orig[2].y-orig[0].y, orig[0].y);
// now one that maps the final optimized control points to the unit right triangle
tform c(cpts[1].x-cpts[0].x, cpts[2].x-cpts[0].x, cpts[0].x,
        cpts[1].y-cpts[0].y, cpts[2].y-cpts[0].y, cpts[0].y);
// now, to get from the original to the final, apply o^-1, then c;
tform oi;
InvertTrans(oi, o);
//tform temp;
MultiplyTrans(t, c, oi);
PrintTransform(of, t);
}

// Improve the correlation, if possible, by tweaking the transform.  pts are the points
// in the original, and pts are their values.  image is a 4Kx4K matrix of doubles, already
// normalized.  dx and dy are the initial estimates of how to map the original points into
// the array image2.
// returns the best correlation obtained.
double ImproveCorrelation(vector<Point> &Plist, vector<double> &spv, vector<double>image2, 
 double dx, double dy, tform &t, FILE *flog)
{
fprintf(of,"Contains %d pixels\n", Plist.size() );
Normalize(spv);
if (dx != BIG)  // if dx == BIG, start improving from the transform we have
    t = tform(1.0, 0.0, dx, 0.0, 1.0, dy);  // otherwise, create a transform with just dx, dy
   
double best_so_far = 0.0;

TryNewOptimizer(Plist, spv, image2, t, flog);

// Now t is the best transform we can find.
fprintf(of,"Best transform is %9.4f %9.4f %10.2f\n                  %9.4f %9.4f %10.2f\n", 
 t.t[0], t.t[1], t.t[2], t.t[3], t.t[4], t.t[5] );
return best_so_far;
}

void WriteTransform(const char *s, tform &t)
{
fprintf(of,"%s %9.4f %9.4f %10.2f\n%s %9.4f %9.4f %10.2f\n", s,
 t.t[0], t.t[1], t.t[2], s, t.t[3], t.t[4], t.t[5] );
}

int iabs(int x){return x >= 0 ? x : -x;}

double Distance(Point &p1, Point &p2)
{
double dx = p1.x - p2.x;
double dy = p1.y - p2.y;
return sqrt(dx*dx + dy*dy);
}

int min(int x, int y){return x < y ? x : y;}
int max(int x, int y){return x > y ? x : y;}


double sgn(double a){return a > 0.0 ? 1.0 :(a < 0.0 ? -1.0 : 0.0);}

// if the original is bigger than 2K, create a downsampled copy (but keep the original)
void Picture::DownsampleIfNeeded()
{
scale = 1;
raster = original;                   // orig always points to a full size raster
if(w <= 2048 && h <= 2048)
    return;

// OK, find a scale factor that will make each dimension <= 2048
int origw = w, origh = h;
do {
    scale *= 2;
    h /= 2;
    w /= 2;
    }
while (w > 2048 || h > 2048);
fprintf(of,"Will scale by %d\n", scale);
if (h*scale != origh || w*scale != origw) {
    fprintf(of,"Image does not divide evenly!\n");
    exit(45);
    }
raster = (uint8 *)malloc(w*h*sizeof(uint8));
int npq = scale*scale;
for(int ix = 0; ix < w; ix++) {
    for(int iy=0; iy < h; iy++) {
	int sum = 0;
        for(int dx=0; dx<scale; dx++) 
	    for(int dy=0; dy<scale; dy++)
		sum += original[ix*scale+dx + (iy*scale+dy)*origw];
        raster[ix+w*iy] = (sum+npq/2)/npq; // rounding
	}
    }
tr.t[2] /= scale;
tr.t[5] /= scale;
}


// return the distance between the point (x2,y) and the line segment from (x0,y0) to (x1,y1)
double LinePointDist(int x0, int y0, int x1, int y1, int x2, int y2)
{
// normalize so (x0.y0) is at (0,0)
x1 -= x0; y1 -= y0;
x2 -= x0; y2 -= y0;
// find parameterized intersection point
double alpha = double(x1*x2 + y1*y2)/double(x1*x1+y1*y1);
// make sure it does nto go off the end of the line
alpha = fmax(alpha,0.0);
alpha = fmin(alpha,1.0);
double dx = x2-alpha*x1;
double dy = y2-alpha*y1;
return sqrt(dx*dx+dy*dy);
}

double LinePointDist(vertex &v0, vertex &v1, vertex &v2)
{return LinePointDist(v0.x, v0.y,  v1.x, v1.y,  v2.x, v2.y);}

// remove from the map all pixels within distance 'dist' of vertex
void RemoveFromMap(vector<uint8> &map, int w, int h, vertex v, double dist)
{
int d = int(dist+1.0);
for(int x = v.x-d; x <= v.x+d; x++) {
    for(int y=v.y-d; y <= v.y+d; y++) {
        if (x < 0 || x >= w || y < 0 || y >= h)
	    continue;
        if ((v.x-x)*(v.x-x) + (v.y-y)*(v.y-y) <= d*d)
	    map[x+y*w] = 0;
	}
    }
}

double DistSquared(vertex &a, vertex &b)
{
double dx = a.x - b.x;
double dy = a.y - b.y;
return dx*dx + dy*dy;
}

// Tells if point c is on the left side of the vector a->b by using the cross product
bool LeftSide(const vertex &a, const vertex &b, const vertex &c)
{
return (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x) > 0;
}

bool LinesCross(const vertex &p1, const vertex &p2, const vertex &p3, const vertex &p4)
{
int a = p2.x - p1.x;
int b = p3.x - p4.x;
int c = p2.y - p1.y;
int d = p3.y - p4.y;
int det = a*d - b*c;

int e = p3.x - p1.x;
int f = p3.y - p1.y;
double alpha = double(d*e - b*f) / double(det);
double beta = double((-c)*e + a*f) / double(det);
bool cross = 0.0 < alpha && alpha < 1.0 && 0.0 < beta && beta < 1.0;
if (cross) {
    fprintf(of,"Aha!  got crossing (%d %d) to (%d %d) crosses (%d %d) to (%d %d), alpha=%f beta=%f\n",
     p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y, alpha, beta);
    double ax = alpha*p2.x + (1-alpha)*p1.x;
    double ay = alpha*p2.y + (1-alpha)*p1.y;
    double bx =  beta*p4.x + (1-beta) *p3.x;
    double by =  beta*p4.y + (1-beta) *p3.y;
    fprintf(of,"Intersect (%f %f) and (%f %f)\n", ax, ay, bx, by);
    }
return cross;
}

// returns TRUE if the line from mid->c crosses any other edges
bool AnyCross(vector<lineseg> &es, vertex mid, vertex c)
{
for(int i=0; i<es.size(); i++) {
    if (LinesCross(es[i].v[0], es[i].v[1], mid, c))
	return true;  // there is a crossing
    }
return false; // there is no crossing
}

// returns the number of crossings of the line from mid->c
int CountCrossings(vector<lineseg> &es, vertex mid, vertex c)
{
int m=0;
for(int i=0; i<es.size(); i++) {
    if (LinesCross(es[i].v[0], es[i].v[1], mid, c))
	m++;
    }
return m;
}

// classes for creating a queue
class qe {
  public:
    qe(int t, double c){to = t; cost = c;}
    bool operator<(const qe &a) const {return cost > a.cost;};  // priority is less if cost is higher
    int to;  // we know how to get to the node 'to' for total cost 'cost'
    double cost;
    };


class gr{
  public:
    gr(int xx, int yy, int b, double c){x = xx; y = yy; back = b; cost = c;}
    int x, y, back;
    double cost;
    };

double LinePointDist(gr &v0, gr &v1, gr &v2)
{return LinePointDist(v0.x, v0.y,  v1.x, v1.y,  v2.x, v2.y);}

int dxs[8] = {1, 1, 0, -1, -1, -1,  0,  1};
int dys[8] = {0, 1, 1,  1,  0, -1, -1, -1};

// try to draw a line from entry s to entry e of vertices.  The cost of this
// approximation is the sum of all the errors at intermediate vertices.
double CostOfApprox(int s, int e, vector<gr> &v)
{
double sum = 0.0;
for(int i=s+1; i != e; i = i+1) 
    sum += LinePointDist(v[s], v[e], v[i]);
return sum;
}

double AreaOfTriangle(const vertex &v0, const vertex &v1, const vertex &v2)
{
return abs((v2.x-v0.x)*(v1.y-v0.y) - (v1.x-v0.x)*(v2.y-v0.y))/2.0;
}
double AreaOfTriangle(const Point &p0, const Point &p1, const Point &p2)
{
return abs((p2.x-p0.x)*(p1.y-p0.y) - (p1.x-p0.x)*(p2.y-p0.y))/2.0;
}

// proposed triangle goes va->vb->vc->va.  See if any vertices are internal to it.
bool AnyInside(vertex &va, vertex &vb, vertex &vc, vector<lineseg> &es, vector<vertex> &ips)
{
for(int i=0; i<es.size(); i++) {
    vertex v = es[i].v[0];
    if (LeftSide(va, vb, v ) && LeftSide(vb, vc, v) && LeftSide(vc, va, v) )
	return true;
    v = es[i].v[1];
    if (LeftSide(va, vb, v ) && LeftSide(vb, vc, v) && LeftSide(vc, va, v) )
	return true;
    }
for(int i=0; i<ips.size(); i++) {
    if (LeftSide(va, vb, ips[i] ) && LeftSide(vb, vc, ips[i]) && LeftSide(vc, va, ips[i]) )
	return true;
    }
return false; // none are inside.
}
void WriteMonochromeTIFF(const char *filename, uint8* buffer, int w, int h, int BitsPerPixel);

// given a vector of integer-valued points, create a bounding polygon that approximates it.
// We can afford a few pixel of errors (perhaps 5 pixels) and do not want too fine of a 
// fragmentation of edges.
// Returns a vector of control points and triangles.  Status = 0 if OK, >=1 if error encountered.
// Also returns a simplified control points list (always 3 points) with 1 triangle for the whole region
int CreateOutline(vector<Point> &pts, vector<vertex> &ControlPoints, vector<triangle> &tris,
 vector<vertex> &SimplePts, vector<triangle> &SimpleTri)
{
// compute min and max.
int xmin = BIG, ymin = BIG, xmax = -BIG, ymax = -BIG;
for(int j=0; j<pts.size(); j++) {
    xmin = min(xmin, int(pts[j].x));
    xmax = max(xmax, int(pts[j].x));
    ymin = min(ymin, int(pts[j].y));
    ymax = max(ymax, int(pts[j].y));
    }
fprintf(of,"region size is [%d %d] in x, [%d %d] in y\n", xmin, xmax, ymin, ymax);
double lmin = minl;    // the normal settings
double lmax = 2*minl;
if (xmax-xmin < minl || ymax-ymin < minl) {  // use different slicing for long skinny pieces
    lmin = fmin(xmax-xmin, ymax-ymin);      // must make lmin short enough to go across the short way
    lmin = fmin(lmin, fmax(xmax-xmin, ymax-ymin)/2.0);  // and less than half the longer side
                                                        // otherwise a point near the middle
                                                        // has no legal distinations at all
    lmax = fmax(xmax-xmin, ymax-ymin)/3.0;
    lmax = fmax(lmax, 2.0*lmin);  // we want at least a 1:2 range from min to max
    fprintf(of," ---Reducing edges from [%f %f] to [%f %f]\n", minl, 2*minl, lmin, lmax);
    }
bool MakeSimple = true;
if (MakeSimple) {
            triangle t;
            t.v[0] = 0; t.v[1]=1; t.v[2]=2;
            SimpleTri.push_back(t);
            int ix1 = int(xmin), iy1 = int(ymin), ix2 = int(xmax), iy2 = int(ymax); 
            if (xmax-xmin > ymax-ymin) { // horizontal
                vertex v0(ix1, iy1, 0), v1(ix2, iy1, 0), v2((ix1+ix2)/2, iy2, 0);
                SimplePts.push_back(v0); SimplePts.push_back(v1); SimplePts.push_back(v2);
		}
            else {    // vertical
                vertex v0(ix1, iy1, 0), v1(ix2, (iy1+iy2)/2, 0), v2(ix1, iy2, 0);
                SimplePts.push_back(v0); SimplePts.push_back(v1); SimplePts.push_back(v2);
		}
            fprintf(of, "Computed simple triangle\n");
    }
// create a bit-map image
int w = xmax-xmin+3;  // a border of blanks all around
int h = ymax-ymin+3;
vector<uint8> map(w*h,0);  // create a map and fill it with zeros
for(int j=0; j<pts.size(); j++) {
    int x = int(pts[j].x) - xmin + 1;
    int y = int(pts[j].y) - ymin + 1;
    map[y*w+x] = 1;
    }
// start from the middle left middle, progress till we hit something
// We work in map coordinates since that is easier
int y0 = h/2;
int x0 = 0;
for(x0=0; x0<w && map[w*y0+x0] == 0; x0++)
    ;
if (x0 == w) {  // this one could possibly happen if the overlap cut the connected region
                // into two parts, neither of which intersects the horizonal midline.
                // Return failure if it does.
    fprintf(of,"Can't find starting point??\n");
#ifdef WRITE_DEBUG_IMAGES
    WriteMonochromeTIFF("Bogus3.tif", &map[0], w, h, 8);
#endif
    return 3;
    //exit(46);
    }
// Now work our way around CCW, till we get back the the beginning
vector<vertex> vertices;
int x = x0, y = y0, dir = 0;
for(; !(x == x0 && y == y0 && dir != 0); ) {
    //fprintf(of,"-- %6d %6d %3d\n", x, y, dir);
    int op = (dir+4)%8;  // opposite to incoming direction
    int k;
    for(k=1; k<=8; k++) {  // one of these directions MUST be non-zero
        int j = (op+k)%8;
        int dx = dxs[j];
	int dy = dys[j];
	if (map[w*(y+dy)+x+dx]) {
            vertex vtx(x,y,j);
            vertices.push_back(vtx);
	    x = x+dx;
	    y = y+dy;
	    dir = j;
	    break;
	    }
	}
    if (k > 8) {
	fprintf(of, "Input was an isolated pixel??\n");
#ifdef WRITE_DEBUG_IMAGES
        WriteMonochromeTIFF("Bogus4.tif", &map[0], w, h, 8);
#endif
	return 4;
	}
    }

// Now compress all the identical edges.  Only save vertices whose direction is different than previous one
vector<vertex> small;
int N = vertices.size();
for(int i=0; i<N; i++) {
    int prev = (i+N-1)%N;  // since using mod N, adding N-1 is same as subtracting 1
    //fprintf(of,"[%d].dir = %d, prev [%d].dir = %d\n", i, vertices[i].dir, prev, vertices[prev].dir);
    if (vertices[i].dir != vertices[prev].dir) {
	vertex v = vertices[i];
        v.orig = i;       // was originally the ith vertex
        small.push_back(v);
	}
    }
fprintf(of,"Had %d original vertices, now %d\n", N, small.size() );
for(int i=0; i<small.size(); i++) {
    int j = small[i].orig;
    fprintf(of," #%3d, was %3d: %8d %6d\n", i, j, small[i].x, small[i].y);
    }

// Now, if any edges are too long, split them.  It's OK to oversplit them since the dynamic
// programming can cope well with this.  (In fact this is helpful, since it gives several
// options through each point.)
for(int i=0; i<small.size(); i++) {
    int next = (i+1)%small.size();
    double dx = small[next].x - small[i].x;
    double dy = small[next].y - small[i].y;
    double len = sqrt(dx*dx + dy*dy);
    //fprintf(of,"len=%f\n", len);
    if (len > lmin/2) {
        int ins = int(len/lmin*2.0)+1; // how many to insert
        for(int j=1; j<=ins; j++) {
	    double frac = double(j)/(ins+1);
            vertex nv(int(small[i].x+dx*frac), int(small[i].y+dy*frac), 0);
            fprintf(of,"Inserting at (%d %d)\n", nv.x, nv.y);
            vector<vertex>::iterator it = small.begin();
	    small.insert(it+i+j, nv);
	    }
	}
    }

for(int i=0; i<small.size(); i++)
    fprintf(of,"--After split edges: vertex %d = (%d %d)\n", i, small[i].x, small[i].y);

// We loop here since we want to use min_l, but sometimes that results in no solutions.
// So we'll try it, and only if it does not work reduce it.
vector<lineseg> edges;
for(int n_iter=1; n_iter <= 10 && edges.size() < 3; n_iter++, lmin = lmin*0.8) {
    fprintf(of,"Pass %d: min length %f, max length %f\n", n_iter, lmin, lmax);
    edges.clear();

    // set all back pointers to -1, and all costs to infinity
    vector<gr> graph;
    for(int i=0; i<small.size(); i++)
	graph.push_back(gr(small[i].x, small[i].y, -1, double(BIG)));

    // now duplicate the first entry after the last.  This will be our target
    graph.push_back(graph[0]);
    // now find the lowest cost path from the first to last vertex.
    // Push the [0] entry on the queue
    priority_queue<qe> q;
    q.push(qe(0, 0.0));
    graph[0].cost = 0.0;

    bool big_print = false;
    while (!q.empty()) {
	qe tp = q.top(); // highest priority == lowest cost
	q.pop();
	fprintf(of,"can get from node 0 to %d for a total cost of %f\n", tp.to, tp.cost);
	if (tp.cost > graph[tp.to].cost) // we already knew how to get here for cheaper
	    continue;
	if (tp.to == graph.size()-1)
	    break;
	// we need to examine all the fanouts of this node
	int s = tp.to;  // the source node
	for(int j=s+1; j<graph.size(); j++) {
	    double dx = graph[j].x - graph[s].x;
	    double dy = graph[j].y - graph[s].y;
	    double len = sqrt(dx*dx + dy*dy);
	    if (big_print) fprintf(of," j=%d len=%f\n", j, len); 
	    if (lmin <= len && len <= lmax) {
		double cost =  graph[s].cost + CostOfApprox(s, j, graph);
		if (big_print) fprintf(of, " src %f, tot %f, tar %f\n", graph[s].cost, cost, graph[j].cost);
		if (cost < graph[j].cost) {// we've found a better way to get there
		    graph[j].cost = cost;
		    graph[j].back = s;
		    if (big_print)fprintf(of," Pushing node %d, cost %f, back %d\n", j, cost, s);
		    q.push(qe(j, cost));
		    }
		}
	    }
	}

    if (graph[graph.size()-1].cost == BIG) {
	fprintf(of,"No path to final vertex??  Will reduce lmin and try again.\n");
	continue;
	}
    //backtrack from the final vertex, adding output edges as we go
    for(int i=graph.size()-1; i != 0; i = graph[i].back) {
	//vertex v(graph[i].x, graph[i].y, 0);
	//RemoveFromMap(map, w, h, v, (lmin+lmax)/2);   // remove all points near existing points
	gr prev = graph[graph[i].back];
	double d=sqrt(pow(double(prev.x-graph[i].x),2.0) + pow(double(prev.y-graph[i].y),2.0) );
	fprintf(of,"edge from (%d %d) to (%d %d),length %7.2f\n",prev.x, prev.y, graph[i].x, graph[i].y, d); 
	edges.push_back(lineseg(prev.x, prev.y, graph[i].x, graph[i].y));
	}
    if (edges.size() <= 2)
	fprintf(of,"Too few edges (%d) for triangulation.  Will reduce lmin and try again.\n", edges.size());
    }

// Now add any internal vertices. Only need one end of each edges since it forms a loop.
for(int i=0; i<edges.size(); i++) {
    RemoveFromMap(map, w, h, edges[i].v[0], (lmin+lmax)/2);
    }
vector<vertex> internal;
for(int i=0; i<w*h; i++) {
     if(map[i]) {
        int y = i/w;
        int x = i-y*w;
        fprintf(of,"Add internal vertex at (%d %d)\n", x, y);
	vertex newv(x,y,0);
        RemoveFromMap(map, w, h, newv, (lmin+lmax)/2);
        vertex outside(-10, y, 0); // vertex has negative x and hence must be outside
        int m = CountCrossings(edges, newv, outside);
        if (m & 1)
            internal.push_back(newv);
        else {
             fprintf(of, "Warning: Supposedly internal point is outside - %d crossings.  Continuing\n", m);
             }
        }
    }

// now divide into triangles.  First, compile a list of (ordered) exterior edges.

// Now, take each edge one at a time.  Find the best triangle for that edge, then remove
// triangle from the figure.  Repeat until no edges left.  Result will be a list of
// control points, and triangles that refer to them.
ControlPoints.clear();
tris.clear();
for(; edges.size() > 0; ) {
    // print all edges in matlab format
    // also, compute the area of the polygon.
    double area = 0.0;
    fprintf(of,"\n--------- Start matlab format ------------------\n");
    for(int j=0; j<edges.size(); j++) {
	fprintf(of,"x=[%d %d]; y=[%d %d]; plot(x,y); hold on;", edges[j].v[0].x, edges[j].v[1].x, edges[j].v[0].y, edges[j].v[1].y);
        if (j > 0 && (j%2) == 0)
	    fprintf(of,"\n");
        area += (edges[j].v[0].x - edges[j].v[1].x)*(edges[j].v[0].y+edges[j].v[1].y)/2.0;
	}
    fprintf(of,"\n");
    fprintf(of," %d edges, %f total area\n", edges.size(), area);
    if (area < 0.0) {
	fprintf(of, "Internal error!  Combined area should not be < 0.\n");
        exit(-1);
	}
    // we'll use the last edge, since it's the easiest one to delete.
    vertex va(edges[edges.size()-1].v[0]);
    vertex vb(edges[edges.size()-1].v[1]);
    edges.pop_back();
    vertex vmid((va.x+vb.x)/2, (va.y+vb.y)/2, 0); // midpoint
    fprintf(of,"working on (%d %d) -> (%d %d)\n", va.x, va.y, vb.x, vb.y);
    // find the closest vertex that's on the right side.
    int id = -2;   // will be -1 if it's an internal vertex, 0 for v[0] of an edge, 1 if v[1] of an edge
    int ii = -1;   // the index into either edges or internal
    double bestd = BIG;
    for(int j=0; j<edges.size(); j++) {
        vertex vc(edges[j].v[0]);
        double d = DistSquared(vmid, vc);
        //fprintf(of,"#1 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if (LeftSide(va, vb, vc) && d < bestd && (!AnyCross(edges, va, vc)) && (!AnyCross(edges, vb, vc))  && (!AnyInside(va, vb, vc, edges, internal)) && AreaOfTriangle(va, vb, vc) > MinTriangleArea) {
	    ii = j; id = 0; bestd = d;
	    }
        vc = edges[j].v[1];
        d = DistSquared(vmid, vc);
        //fprintf(of,"#2 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if (LeftSide(va, vb, vc) && d < bestd && (!AnyCross(edges, va, vc)) && (!AnyCross(edges, vb, vc))  && (!AnyInside(va, vb, vc, edges, internal)) && AreaOfTriangle(va, vb, vc) > MinTriangleArea) {
	    ii = j; id = 1; bestd = d;
	    }
        }
    for(int j=0; j<internal.size(); j++) {
        vertex vc(internal[j]);
        double d = DistSquared(vmid, vc);
        //fprintf(of,"#3 (%d %d) dist %f %d\n", vc.x, vc.y, d, LeftSide(va, vb, vc) );
        if (LeftSide(va, vb, vc) && d < bestd && (!AnyCross(edges, va, vc)) && (!AnyCross(edges, vb, vc))  && (!AnyInside(va, vb, vc, edges, internal)) && AreaOfTriangle(va, vb, vc) > MinTriangleArea) {
	    ii = j; id = -1; bestd = d;
	    }
        }
    if (id == -2) {
	fprintf(of, "No legal triangle at all?? %d triangles so far, area limit %.1f.\n", tris.size(), MinTriangleArea );
	if (tris.size() > 0)
	    break;     // at least some triangles generated.  Stop here and say we are OK.
	else { // none at all.  Make one up from the bounding box
            triangle t;
            t.v[0] = 0; t.v[1]=1; t.v[2]=2;
            tris.push_back(t);
            int ix1 = int(xmin), iy1 = int(ymin), ix2 = int(xmax), iy2 = int(ymax); 
            if (xmax-xmin > ymax-ymin) { // horizontal
                vertex v0(ix1, iy1, 0), v1(ix2, iy1, 0), v2((ix1+ix2)/2, iy2, 0);
                ControlPoints.push_back(v0); ControlPoints.push_back(v1); ControlPoints.push_back(v2);
		}
            else {    // vertical
                vertex v0(ix1, iy1, 0), v1(ix2, (iy1+iy2)/2, 0), v2(ix1, iy2, 0);
                ControlPoints.push_back(v0); ControlPoints.push_back(v1); ControlPoints.push_back(v2);
		}
            fprintf(of, "STAT: fallback triangle from %d points\n", pts.size() );
	    return 0;  // 1 legal triangle
            }
	}
    vertex vc = id < 0 ? internal[ii] : edges[ii].v[id];
    // now the triangle goes va->vb->vc->va
    fprintf(of,"Triangle (%d %d) (%d %d) (%d %d) area %f\n", va.x, va.y, vb.x, vb.y, vc.x, vc.y, AreaOfTriangle(va, vb, vc));
    triangle t;
    vector<vertex> temp; temp.push_back(va); temp.push_back(vb); temp.push_back(vc);
    for(int j=0; j<3; j++) {
	int k;
	for(k=0; k<ControlPoints.size(); k++)
	    if (ControlPoints[k] == temp[j]) {
		t.v[j] = k;
		break;
	    }
        if (k == ControlPoints.size()) { // did not find it; add it
            ControlPoints.push_back(temp[j]);
            t.v[j] = k;
	    }
        }
    tris.push_back(t);
    // for the edges vb->vc and vc->va, if they are in the list of edges, remove them
    // if they are not in the list, then add them.
    lineseg newone(vb,vc);
    int start = edges.size();
    for(int j=0; j<start; j++) {
	if (edges[j] == newone) {
	    edges.erase(edges.begin()+j);
	    break;
	    }
        }
    if(edges.size() == start)
	edges.push_back(lineseg(vc,vb));
    lineseg newtwo(vc,va);
    start = edges.size();
    for(int j=0; j<start; j++) {
	if (edges[j] == newtwo) {
	    edges.erase(edges.begin()+j);
	    break;
	    }
        }
    if(edges.size() == start)
	edges.push_back(lineseg(va,vc));
    // if the edge came from the internal list, it's internal no longer, so remove it.
    if (id == -1)
	internal.erase(internal.begin() + ii);
    }
fprintf(of,"STAT: from %d pts, got %d triangles using %d control points\n", pts.size(), tris.size(), ControlPoints.size() );

// Convert back to input coordinates
for(int j=0; j<ControlPoints.size(); j++) {
    ControlPoints[j].x += (xmin - 1 );
    ControlPoints[j].y += (ymin - 1 );
    }
return 0;
}

// Find the best triangle.  If we are inside it's obvious; otherwise pick the closest
int BestTriangle(Point &pt, vector<triangle> &tris, vector<vertex> &ControlPoints)
{
vertex v(int(pt.x), int(pt.y), 0);
for(int i=0; i<tris.size(); i++) {
    vertex v0 = ControlPoints[tris[i].v[0]], v1 = ControlPoints[tris[i].v[1]], v2 = ControlPoints[tris[i].v[2]];
    if (LeftSide(v0, v1, v) && LeftSide(v1, v2, v) && LeftSide(v2, v0, v) )
	return i;
     }
//OK, it's not inside any of them.  Find the closest one...
double dmin = BIG;
int best = -1;
for(int i=0; i<tris.size(); i++) {
    vertex v0 = ControlPoints[tris[i].v[0]], v1 = ControlPoints[tris[i].v[1]], v2 = ControlPoints[tris[i].v[2]];
    double d = LinePointDist(v0, v1, v);
    if (d < dmin)
	{dmin = d; best = i;}
    d = LinePointDist(v1, v2, v);
    if (d < dmin)
	{dmin = d; best = i;}
    d = LinePointDist(v2, v0, v);
    if (d < dmin)
	{dmin = d; best = i;}
    }
return best;
}


void FindConnRegions(vector<ConnRegion> &cr, uint8* fold_mask, int w, int h, int scale)
{
fprintf(of,"Finding connected regions: w=%d h=%d scale=%d\n", w, h, scale);
int n = w * h;
int big = 0;
for(int i=0; i<n; i++)
    big = max(big, fold_mask[i]);
cr.clear();
cr.insert(cr.begin(), big, ConnRegion());
for(int i=1; i <= big; i++) {
    cr[i-1].id = i;
    vector<vertex> vv;
    for(int j=0; j<n; j++) {
	if (fold_mask[j] == i) {
	    int iy = j/w;
            int ix = j - iy*w;
            vv.push_back(vertex(ix/scale, iy/scale, 0));
	    }
       }
    // may contain duplicates from down-scaling.  Sort, then take unique values.
    sort(vv.begin(), vv.end());
    cr[i-1].xmin = cr[i-1].ymin = BIG;
    cr[i-1].xmax = cr[i-1].ymax = -BIG;
    for(int j=0; j<vv.size(); j++) {
	if (j == 0 || !(vv[j] == vv[j-1])) {
	    cr[i-1].pts.push_back(Point(vv[j].x, vv[j].y));
            cr[i-1].xmin = min(cr[i-1].xmin, vv[j].x);
            cr[i-1].xmax = max(cr[i-1].xmax, vv[j].x);
            cr[i-1].ymin = min(cr[i-1].ymin, vv[j].y);
            cr[i-1].ymax = max(cr[i-1].ymax, vv[j].y);
	    }
	}
    fprintf(of, "CR %d, %d points, (%d %d) to (%d %d)\n", i, cr[i-1].pts.size(), cr[i-1].xmin, cr[i-1].ymin, cr[i-1].xmax, cr[i-1].ymax);
    }
// drop any that do not have enough points
for(int i=0; i<cr.size(); i++) {
    if (cr[i].pts.size() < 0.9 * MinMapArea) {
	cr.erase(cr.begin()+i); // delete the entry
	i--;  // will need to look at this index again
	}
    }
fprintf(of, "%d CRs were big enough\n", cr.size() );
}

bool BigEnough(int sx, int sy, void *a)
{
long min_area = long(a);
return sx*sy > min_area;
}
bool EnoughPoints(int count1, int count2, void *a)
{
long min_count = long(a);
//printf("c1, c2, m = %10d %10d %10ld\n", count1, count2, min_count);
return count1 > min_count && count2 > min_count;
}

// replaces vectors of Points and values with decimated vectors.
void DecimateVector(vector<Point> &p, vector<double> &v, int n)
{
vector<Decimate> dv;
for(int i=0; i<p.size(); i++) {
    Decimate d;
    d.x = int(floor(p[i].x/n));
    d.y = int(floor(p[i].y/n));
    d.v = v[i];
    dv.push_back(d);
    }
// now sort it
sort(dv.begin(), dv.end());
//Now create new lists
p.clear();
v.clear();
int j;
for(int i=0; i<dv.size(); i = j) {
    double sum = dv[i].v;
    for(j=i+1; j < dv.size() && !(dv[i] < dv[j]); j++)
	sum += dv[j].v;
    double avg = sum/(j-i);
    p.push_back(Point(dv[i].x, dv[i].y));
    v.push_back(avg);
    }
}

void JustTesting(vector<Point> ap, vector<double> av, vector<Point> bp, vector<double> bv, FILE *of)
{
long size = 100000;
int mult = 1;
vector<CD> ftc;
for(int i=0; i<8; i++) {
    double dx, dy;
    ftc.clear(); // no caching wanted here
    double c = FindNormCorrelation(ap, av, bp, bv, 
     dx, dy, 0, 0, 4000, of, BigEnough, (void *)size, NULL, NULL, ftc);
    fprintf(of, " JT: %4d by %4d, size %6d, corr %f, dx, dy = %9.2f, %9.2f\n", 
     int(sqrt(ap.size())), int(sqrt(bp.size())), 
     size, c, dx*mult, dy*mult);
    DecimateVector(ap, av, 2);
    DecimateVector(bp, bv, 2);
    size = size/4;
    mult = mult * 2;
    }
}
// Tries to find a mapping from one connected region to another.
// a = above, b = below, we look for a mapping from a->b
// The picture has the whole image, the connected region tells what part we are concerned with.
// Returns TRUE and adds one or more entries to GUESS if it finds a match.  Otherwise returns false
bool ApproximateMatch(Picture &pa, ConnRegion &a, Picture &pb, ConnRegion &b, vector<tform> &guesses,
FILE *flog)
{

// now we need to find how each of these regions map to the picture below
// if indeed they do at all.
int w = pa.w;		// must be the same for a and b, verified earlier
int h = pa.h;
int sc = pa.scale;
int fullw = pa.w*sc;
int fullh = pa.h*sc;
int npixels = b.pts.size();
fprintf(of," Looking for an approximate match - scale = %d\n", sc);
MeanStd mm;
for(int i=0; i<npixels; i++) {
    int ix = int(b.pts[i].x);
    int iy = int(b.pts[i].y);
    mm.Element(pb.raster[iy*w+ix]);
    }
fprintf(of,"below: %d real pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/(w*h));
double mean2, std2;
mm.Stats(mean2, std2);  // find existing statistics
fprintf(of,"On the target image,  mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy, with image in lower left 2Kx2K, so FFT will work.  
// Background pixels map to zero.
// Also, create a map that just tells whether each point is in the region, or not
vector<double> bv;
for(int i=0; i<npixels; i++) {
    int x = int(b.pts[i].x);
    int y = int(b.pts[i].y);
    double pix = (pb.raster[y*w + x] - mean2)/std2;
    //if (abs(pix) > 2.5)    // outlying pixel value
	//pix = 0.0;  // will be set to the mean value (0)
    bv.push_back(pix);
    }

// Now try a normalized correlation to make a quick estimate
vector<double> av;
for(int k=0; k<a.pts.size(); k++) {
    int x = int(a.pts[k].x);
    int y = int(a.pts[k].y);
    av.push_back(pa.raster[y*w+x]);
    }

// now create a smaller set of images.  We want #pixels to be 40,000 < m < 160,000, by experiment
vector<Point> apts = a.pts;
vector<Point> bpts = b.pts;
int mult = 1;
long npatch = int(NumCorrPatchSize);  // number of pixels for a match
//for(; apts.size() > 160000 && bpts.size() > 160000; ) {
for(int j=0; j<3; j++) {  //always use 3 rounds of reduction, for now
    DecimateVector(apts, av, 2);
    DecimateVector(bpts, bv, 2);
    npatch = npatch/4;
    mult = mult * 2;
    }
fprintf(of, "After reduction %d %d pts, patch %d, mult %d\n", apts.size(), bpts.size(), npatch, mult);

Normalize(av);  // normalize, to make it 0 mean, so the padding on the rotated rectangles
                // will not have too big of an effect

// Here we try a bunch of angles to try to find a match on the reduced size image

double span = 90.0;  // total span searched, in degrees
const double RES  = 0.01;  // angular resolution

int ASIZE = int(span/180.0*PI/RES+0.5); // span in radians
ASIZE = (ASIZE+1)/2; ASIZE *= 2;    // make sure it is even
int MID = ASIZE/2;                  // (0..MID-1) center (MID..ASIZE-1)

vector<double> angs(ASIZE, 0.0);
vector<Point>    ds(ASIZE, Point(0.0,0.0));  // deltas in x,y for each match
vector<tform>   tfs(ASIZE);                  // transform used for each match (may include scaling, skew, etc.
vector<double> corr(ASIZE, -10.0);  /// -10 is an impossible value for a normalized correlation
// set the angles.  We do it this way so no angle is exactly 'center', since that gives anamolous results
// since it does not interpolate if it exactly 0,90,180,270(unlike all other angles).
for(int i=0; i<MID; i++) {
    angs[MID+i] =   AngleCenter/180.0*PI + RES/2 + i*RES;
    angs[MID-1-i] = AngleCenter/180.0*PI - RES/2 - i*RES;
    }

bool found;
vector<int> great;  // at the end of the loop, the high value
vector<CD> ftc;  // fourier transform cache
int next = 1;
double scale = ApproxScale;
tform skew(1.0, 0.0, 0.0, InitialSkew, 1.0, 0.0);

// Loop over all angles, starting at the midde and working out.  Stop when we reach a 'great'
// correlation, or when the array is full.
for(; great.size() == 0 && next != -1;){
    // we want an entry that is over theshold, with two neighbors each side defined, and bigger than each
    great.clear();
    next = -1; // next one to be evaluated.
    // this loop has two functions.  Look for 'great' correlations, and look for a next position to be
    // tried.  We want one next to a 'good' correlation, if we can find one
    for(int i=2; i<ASIZE-2 && great.size() == 0; i++) {
        if(corr[i] >= MinNormCorrelation && corr[i] >= corr[i-1] && corr[i] >= corr[i+1] &&
         corr[i+1] >= corr[i+2] && corr[i-1] >= corr[i-2]) {
	    if (corr[i-1] < -2.0)
		next = i-1;
            else if (corr[i+1] < -2.0)
	        next = i+1;
            else if (corr[i-2] < -2.0)
	        next = i-2;
            else if (corr[i+2] < -2.0)
	        next = i+2;
            else if (corr[i] > MinNormGreat)
	        great.push_back(i);
	    }
        }
    // if 'next' has not been set yet, pick one.  Start from the center and work out.
    if (next == -1 && great.size() == 0) {
        for(int i=0; i < MID; i++) {
            if (corr[MID-1-i] < -2.0)
                {next = MID-1-i; break;}
            if (corr[MID+i] < -2.0)
                {next = MID+i; break;}
	    }
        }
    if (great.size() == 0 && next != -1) {
	// evaluate entry 'next'
	double theta = angs[next];
        printf(" Next = %d, ang=%f\n", next, theta);
	vector<Point> ps;
	tform ao; // angle only
	ao.t[0] = InitXScale*scale*cos(theta); ao.t[1] = -InitYScale*scale*sin(theta);
	ao.t[3] = InitXScale*scale*sin(theta); ao.t[4] =  InitYScale*scale*cos(theta);
        tform use;
        MultiplyTrans(use, ao, skew);
	for(int k=0; k<apts.size(); k++) {
	    Point p = apts[k];
	    Transform(p, use);
	    ps.push_back(p);
	    }
	corr[next] = FindNormCorrelation(ps, av, bpts, bv, ds[next].x, ds[next].y, 0, 0, 4000, of, BigEnough, (void *)npatch, EnoughPoints, (void *)npatch, ftc);
        tfs[next] = use; tfs[next].t[2] = ds[next].x; tfs[next].t[5] = ds[next].y;  // save the transform of best match
	fprintf(of,"New %.5f at angle %8.3f---> %8.2f %8.2f %c\n", corr[next], theta, ds[next].x, ds[next].y, 
         corr[next] > MinNormCorrelation ? '*' : ' ');
        }
    }
vector<int> good;
double limit1 = TransAreClose/mult; //the translations at neighboring angles a[i-1] and a[i+1] should be this close
double limit2 = 2.0*limit1;         //same for those that are 2 away
if (great.size() == 0 ) {  // no great candidates.  Look for good ones.
    MeanStd mm;
    for(int i=2; i<ASIZE-2; i++) 
	mm.Element(corr[i]);
    double mean, std;
    mm.Stats(mean, std);
    if (std < 0.001)
	return false;  // Not enough difference between solutions (most likely all were 0 since
			// no region was big enough
    for(int i=2; i<ASIZE-2 && great.size() == 0; i++) {
        double nval = (corr[i]-mean) / std;
        // a point will be a candidate if (a) it's more than 2 std dev above the mean
        // and (b) there is no nearby point with a higher correlation and a similar offset
        char show;  // will show if it's a candidate
        if (nval < 2.0 ||    // BLOOT - should be programmible
          (Distance(ds[i], ds[i-1]) < limit1 && corr[i-1] > corr[i]) ||
          (Distance(ds[i], ds[i+1]) < limit1 && corr[i+1] > corr[i]) ||
          (Distance(ds[i], ds[i-2]) < limit2 && corr[i-2] > corr[i]) ||
          (Distance(ds[i], ds[i+2]) < limit2 && corr[i+2] > corr[i]) )
	    show = ' ';
	else {
	    good.push_back(i);
	    show = '*';
	    }
	fprintf(of, "In order, ang = %8.4f, dx,dy= %7.1f %7.1f, corr = %7.4f, norm val = %7.3f %c\n", 
         angs[i], ds[i].x, ds[i].y, corr[i], nval, show);
	}
    if (good.size() == 0)  // none of these either
	return false;
    }

vector<int> returns = great.size() > 0 ? great : good;
fprintf(of, "STAT:Returning %d potential solutions\n", returns.size() );
fprintf(flog, "STAT:Returning %d potential solutions\n", returns.size() );
for(int k=0; k<returns.size(); k++) {
    int i = returns[k];
    fprintf(of, "Index %d, angle %f, corr %f \n", i, angs[i], corr[i] );
    double cbest = corr[i];
    double deltax = 0.0;  // start by assuming no interpolation
    // now interpolate to find the best theta, if neighbors have similar offsets
    if (Distance(ds[i], ds[i-1]) <= limit1 && Distance(ds[i], ds[i+1]) <= limit1 ) {
	double mi = corr[i-1] - corr[i];  // relative value at minus 1 coordinate
	double pl = corr[i+1] - corr[i];  // same for +1 coordinate
	double aa = (mi+pl)/2;
	double bb = (pl-mi)/2;  // local equation now a*x^2 + b*x + c (but c == 0) by construction
	deltax = -bb/(2*aa);
	fprintf(of,"Interp: (%f %f %f) %f %f -> dx = %f\n", mi, corr[i], pl, aa, bb, deltax);
	if (deltax < 0) {
	    deltax += 1.0;
	    i--;
	    }
        }
    double bestTheta = angs[i] + deltax*RES;
    fprintf(of,"STAT: case %d, Best correlation is %f, at angle %f \n", k, cbest, bestTheta);
    fprintf(flog,"STAT: case %d, Best correlation is %f, at angle %f \n", k, cbest, bestTheta);
    // now interpolate everything
    tform tr_guess;
    for(int j=0; j<6; j++)
	tr_guess.t[j] = (1-deltax)*tfs[i].t[j] + deltax*tfs[i+1].t[j];
    //tr_guess.t[0] =  InitXScale*scale*((1-deltax)*cos(angs[i]) + deltax*cos(angs[i+1]));
    //tr_guess.t[1] = -InitYScale*scale*((1-deltax)*sin(angs[i]) - deltax*sin(angs[i+1]));
    //tr_guess.t[2] =  (1-deltax)*ds[i].x      + deltax*ds[i+1].x;
    //tr_guess.t[3] =  InitXScale*scale*((1-deltax)*sin(angs[i]) + deltax*sin(angs[i+1]));
    //tr_guess.t[4] =  InitYScale*scale*((1-deltax)*cos(angs[i]) + deltax*cos(angs[i+1]));
    //tr_guess.t[5] =  (1-deltax)*ds[i].y      + deltax*ds[i+1].y;
    tr_guess.t[2] *= mult;
    tr_guess.t[5] *= mult;
    fprintf(of,"Found approximate match; transform is:"); PrintTransform(of, tr_guess);
    // Now try to improve the guess, if possible.
    vector<tform> tweaks(9);
    // tweaks[0] will be the identity matrix
    tweaks[1] = tform(0.995,  0.0,   0.0,  0.0,   0.995, 0.0);  // slightly smaller
    tweaks[2] = tform(1.005,  0.0,   0.0,  0.0,   1.005, 0.0);  // slightly bigger
    tweaks[3] = tform(0.995,  0.0,   0.0,  0.0,   1.005, 0.0);  // X smaller, Y bigger
    tweaks[4] = tform(1.005,  0.0,   0.0,  0.0,   0.995, 0.0);  // slightly bigger
    tweaks[5] = tform(1.000,  0.005, 0.0,  0.005, 1.001, 0.0);  // toe in
    tweaks[6] = tform(1.000, -0.005, 0.0, -0.005, 1.000, 0.0);  // toe out
    tweaks[7] = tform(1.000,  0.005, 0.0, -0.005, 1.001, 0.0);  // Small rotate
    tweaks[8] = tform(1.000, -0.005, 0.0,  0.005, 1.000, 0.0);  // other way
    int works = true;
    for(; works; ) {
        // try all 9 tweaks
        vector<tform> rtran(9);  // resulting transforms
        vector<double> rslts(9);
        int best = -1;
        double best_val = -10.0;
        for(int j=0; j<9; j++) {
	    tform current = tr_guess;
	    current.t[2] = current.t[5] = 0.0;  // remove offset
	    MultiplyTrans(rtran[j], current, tweaks[j]);
	    vector<Point> ps;
	    for(int k=0; k<apts.size(); k++) {
		Point p = apts[k];
		Transform(p, rtran[j]);
		ps.push_back(p);
		}
            // find the offset of the current best match.  We will only look near this spot
	    int ix = int(tr_guess.t[2]/mult);
	    int iy = int(tr_guess.t[5]/mult);
            printf("Initial guess x,y= %d %d\n", ix, iy);
            const int radius = 4000;
            double co, dx, dy;  // correlation and offset
	    co = FindNormCorrelation(ps, av, bpts, bv, dx, dy, ix, iy, radius, of, BigEnough, (void *)npatch, EnoughPoints, (void *)npatch, ftc);
            fprintf(of, "Tweak %d, rslt %f, dx,dy= %f %f\n", j, co, dx, dy);
            rslts[j] = co;
            rtran[j].t[2] = dx*mult;
            rtran[j].t[5] = dy*mult;
            if (co > best_val) {
                best_val = co;
                best = j;
                }
            }
        // now, see who was the best.
        works = best != 0;  // we found one better than identity transform
        if (works) {
	    tr_guess = rtran[best];
            fprintf(of, "New best transform "); PrintTransform(of, tr_guess);
            }
        }
    guesses.push_back(tr_guess);
    }
return true;
}

// Compute a Difference-of-Gaussian raster, if it is not already there.
void Picture::MakeDoGExist(double r1, double r2, vector<CD> &cache)
{
if (DoG.size() > 0)
    return;
DoG.resize(w*h);
bool use_DoG = OuterR > 0.0;  // Makes no sense unless outer R is fairly big
if(!use_DoG) {                // if not using a filter, just copy the array
    for(int k=0; k<w*h; k++)
	DoG[k] = raster[k];
    return;
    }
int N = 1; // size of FFT
while (N < max(w,h)+3*r2)
    N *= 2;
printf("DoG using size %d\n", N);
int M = N*(N/2+1);  // size of array of complex coefficients

if (cache.size() == 0) {
    vector<double> image(N*N, 0.0);
    double sum1 = 0.0, sum2 = 0.0;  // for normalizing.  Though each integral should = 2*pi
       // they are not quite the same due to quantization and removal of tails.
    for(int x = int(-3.0*r2); x <= int(3.0*r2); x++) {
        for(int y = int(-3.0*r2); y <= int(3.0*r2); y++) {
            double rad2 = x*x+y*y;  // radius squared
            sum1 += exp(-(rad2/(2*r1*r1)));
            sum2 += exp(-(rad2/(2*r2*r2)));
	    }
	}
    printf("%f should equal %f\n", sum1/r1/r1, sum2/r2/r2);
    for(int x = int(-3.0*r2); x <= int(3.0*r2); x++) {
        for(int y = int(-3.0*r2); y <= int(3.0*r2); y++) {
            double rad2 = x*x+y*y;
            double v1 = exp(-(rad2/(2*r1*r1)));
            double v2 = exp(-(rad2/(2*r2*r2)));
            int ix = (x >= 0) ? x : N+x;
            int iy = (y >= 0) ? y : N+y;
            image[iy*N+ix] = v1/sum1 - v2/sum2;
	    }
	}
    // Now fft it, and put the result in the cache.
    cache.resize(M);
    fftw_plan p;
    p = fftw_plan_dft_r2c_2d(N, N, (double *)(&image[0]), (double (*)[2])(&cache[0]), FFTW_ESTIMATE);
    fftw_execute(p); // execute the plan
    fftw_destroy_plan(p);
    }

// Now compute the filtered pass image.  First create a copy with a mean of 0, 
// to minimize edge problems
vector<double> t(w*h);
for(int i=0; i<w*h; i++)
    t[i] = raster[i];
Normalize(t);

// Now copy it into a bigger raster
vector<double> imag(N*N, 0.0);
for(int i=0; i<w*h; i++) {
    int iy = i/w;
    int ix = i - iy*w;
    imag[ix + N*iy] = t[i];
    }
vector<CD> freq(M);
fftw_plan p;
p = fftw_plan_dft_r2c_2d(N, N, (double *)(&imag[0]), (double (*)[2])(&freq[0]), FFTW_ESTIMATE);
fftw_execute(p); // execute the plan
fftw_destroy_plan(p);

// convolve my multiplying in the frequency domain
for(int i=0; i<M; i++)
    freq[i] = freq[i] * cache[i];

// Now, fft back
p = fftw_plan_dft_c2r_2d(N, N, (double (*)[2])(&freq[0]), (double *)(&imag[0]), FFTW_ESTIMATE);
fftw_execute(p); // execute the plan
fftw_destroy_plan(p);

vector<double> v(w*h);
for(int i=0; i<w*h; i++) {
    int iy = i/w;
    int ix = i - w*iy;
    v[i] = imag[ix + N*iy];
    }
Normalize(v);
DoG.resize(w*h);
for(int i=0; i<w*h; i++) {
    int pix = int(127 + 32*v[i]);
    pix = min(pix, 255);
    pix = max(pix, 0);
    DoG[i] = pix;
    }
}

double ImproveMesh(vector<Point> &newpts, vector<double> &spv, // describes source image
 vector<vertex> &ControlPoints, vector<triangle> &tris,        // describes triangles in source
 vector<double> &image2,                                       // describes the target image
 tform &tr_guess,                                              // approximate transform
 double threshold,					       // how good results must be before acceptance
 vector<tform> &transforms, vector<Point> &centers,            // results go here
 FILE *flog,
 const char *describe)                                         // string describing operation, for log file
{

// for each triangle, create a matrix that transforms to barycentric coordinates
for(int k=0; k<tris.size(); k++) {
    vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
    fprintf(of,"Tri: (%d %d) (%d %d) (%d %d)\n", v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
    double a[3][3];
    a[0][0] = v0.x; a[0][1] = v1.x; a[0][2] = v2.x;
    a[1][0] = v0.y; a[1][1] = v1.y; a[1][2] = v2.y;
    a[2][0] = 1.0;       a[2][1] = 1.0;       a[2][2] = 1.0;
    Invert3x3matrix(tris[k].a, a);
    }
for(int k=0; k<tris.size(); k++) {  // write in matlab form
    vertex v0 = ControlPoints[tris[k].v[0]], v1 = ControlPoints[tris[k].v[1]], v2 = ControlPoints[tris[k].v[2]];
    fprintf(of,"x=[%d %d %d %d]; y=[%d %d %d %d]; plot(x,y); hold on;\n", v0.x, v1.x, v2.x, v0.x,
     v0.y, v1.y, v2.y, v0.y);
    }
vector<ConnRegion> subs(tris.size());

for(int k=0; k<subs.size(); k++) {
    subs[k].xmin = subs[k].ymin = BIG;
    subs[k].xmax = subs[k].ymax = -BIG;
    }
// assign each point to one of the sub-regions
vector<vector<double> >lambda;
for(int k=0; k<newpts.size(); k++) {
    Point pt = newpts[k];
    int reg = BestTriangle(pt, tris, ControlPoints);
    subs[reg].pts.push_back(pt);
    subs[reg].xmin = min(subs[reg].xmin, int(pt.x));
    subs[reg].ymin = min(subs[reg].ymin, int(pt.y));
    subs[reg].xmax = max(subs[reg].xmax, int(pt.x));
    subs[reg].ymax = max(subs[reg].ymax, int(pt.y));
    // Map the point into barycentric coordinates for its triangle
    //vertex v0 = ControlPoints[tris[reg].v[0]], v1 = ControlPoints[tris[reg].v[1]], v2 = ControlPoints[tris[reg].v[2]];
    //fprintf(of," pt (%f %f) in tri: (%d %d) (%d %d) (%d %d)\n", pt.x, pt.y, v0.x, v0.y, v1.x, v1.y, v2.x, v2.y);
    double l[3];
    l[0] = tris[reg].a[0][0]*pt.x + tris[reg].a[0][1]*pt.y + tris[reg].a[0][2]*1.0;
    l[1] = tris[reg].a[1][0]*pt.x + tris[reg].a[1][1]*pt.y + tris[reg].a[1][2]*1.0;
    l[2] = tris[reg].a[2][0]*pt.x + tris[reg].a[2][1]*pt.y + tris[reg].a[2][2]*1.0;
    //fprintf(of,"prop %f %f %f\n", l[0], l[1], l[2]);
    vector<double>lam(ControlPoints.size(), 0.0);
    for(int m=0; m<3; m++)
	lam[tris[reg].v[m]] = l[m];
    //for(int m=0; m<lam.size(); m++)
	//fprintf(of,"%.4f ", lam[m]);
    //fprintf(of,"\n");
    lambda.push_back(lam);
    }
// optimize the mesh.
// next, transform the control points to the target frame
vector<Point> cpts; // make a real-valued copy
for(int k=0; k<ControlPoints.size(); k++)
    cpts.push_back(Point(ControlPoints[k].x, ControlPoints[k].y));
vector<Point> orig = cpts;
Transform(cpts, tr_guess);
vector<Point> pre = cpts;
// call the optimizer to find a better spot for the control points
double corr = ImproveControlPts(cpts, lambda, spv, image2, 4096, flog, describe, threshold);
if (corr < threshold) {
    return corr;
    }
// write out the original, and the modified, in MatLab format:
for(int k=0; k<tris.size(); k++) {  // write in matlab form
    Point v0 = orig[tris[k].v[0]], v1 = orig[tris[k].v[1]], v2 = orig[tris[k].v[2]];
    fprintf(of,"x=[%f %f %f %f]; y=[%f %f %f %f]; plot(x,y); hold on;\n", v0.x, v1.x, v2.x, v0.x,
     v0.y, v1.y, v2.y, v0.y);
    }
for(int k=0; k<tris.size(); k++) {  // write in matlab form
    Point v0 = cpts[tris[k].v[0]], v1 = cpts[tris[k].v[1]], v2 = cpts[tris[k].v[2]];
    fprintf(of,"x=[%f %f %f %f]; y=[%f %f %f %f]; plot(x,y); hold on;\n", v0.x, v1.x, v2.x, v0.x,
     v0.y, v1.y, v2.y, v0.y);
    }
// see how far the points moved.
for(int k=0; k<cpts.size(); k++) {
    double d = Distance(cpts[k], pre[k]);
    fprintf(of,"%6d: dx=%8.2f dy=%8.2f d=%8.2f\n", k, cpts[k].x-pre[k].x, cpts[k].y-pre[k].y, d);
    }
// check the change in triangle area
double max_change = 0.0;
double sum_old = 0.0, sum_now = 0.0;
for(int k=0; k<tris.size(); k++) {  // write in matlab form
    Point v0 = orig[tris[k].v[0]], v1 = orig[tris[k].v[1]], v2 = orig[tris[k].v[2]];
    double old = AreaOfTriangle(v0, v1, v2);
    if (ApproxScale != 1.0 || InitXScale != 1.0 || InitYScale != 1.0) {
        double f = ApproxScale*ApproxScale*InitXScale*InitYScale;
        fprintf(of, " Modifying old area from %f to %f because scale change was specified\n", old, old*f);
        old *= f;
        }
    sum_old += old;
    v0 = cpts[tris[k].v[0]]; v1 = cpts[tris[k].v[1]]; v2 = cpts[tris[k].v[2]];
    double now = AreaOfTriangle(v0, v1, v2);
    sum_now += now;
    double pct = (now-old)/old*100.0;
    fprintf(of, " Triangle %d, area was %10.1f, is %10.1f, %6.1f%%\n", k, old, now, pct);
    max_change = fmax(max_change, fabs(pct));
    }
double pct = (sum_now-sum_old)/sum_old*100.0;
fprintf(of, " Combined: area was %10.1f, is %10.1f, %6.1f%%\n", sum_old, sum_now, pct);

// see if the area change is reasonable.  If the area is big (>500K pixels) allow 12%
// Otherwise allow 9%.  BLOOT - Should make these parameters
double delta_threshold = sum_now > 500000 ? 12.0 : 9.0;
if (max_change >= 30.0 || fabs(pct) > delta_threshold) {
    fprintf(  of, "STAT: Area change too big (%8.2f%% %8.2f%%), limit %f\n", max_change, pct, delta_threshold);
    fprintf(flog, "STAT: Area change too big (%8.2f%% %8.2f%%), limit %f\n", max_change, pct, delta_threshold);
    return 0.0;
    }
// compute a transform, and a center, for each triangle
for(int k=0; k<tris.size(); k++) {
    int i0 = tris[k].v[0];
    int i1 = tris[k].v[1];
    int i2 = tris[k].v[2];
    // Now, find a transformation that maps ORIG into the new cpts
    // first, create a transform that maps a unit right triangle to the original pts
    tform o(orig[i1].x-orig[i0].x, orig[i2].x-orig[i0].x, orig[i0].x,
     orig[i1].y-orig[i0].y, orig[i2].y-orig[i0].y, orig[i0].y);
    // now one that maps the final optimized control points to the unit right triangle
    tform c(cpts[i1].x-cpts[i0].x, cpts[i2].x-cpts[i0].x, cpts[i0].x,
     cpts[i1].y-cpts[i0].y, cpts[i2].y-cpts[i0].y, cpts[i0].y);
    // now, to get from the original to the final, apply o^-1, then c;
    tform oi,t;
    InvertTrans(oi, o);
    //tform temp;
    MultiplyTrans(t, c, oi);
    PrintTransform(of, t);
    if (abs(t.t[0]-1) > 0.1 || abs(t.t[4]-1) > 0.1) {
	fprintf(of, "Very suspicious transform: %d %d %d\n",
	 tris[k].v[0], tris[k].v[1], tris[k].v[2] );
	fprintf(of, "(%f %f) (%f %f) (%f %f)\n",
	 orig[i0].x, orig[i0].y, orig[i1].x, orig[i1].y,orig[i2].x, orig[i2].y);
	fprintf(of, "(%f %f) (%f %f) (%f %f)\n",
	 cpts[i0].x, cpts[i0].y, cpts[i1].x, cpts[i1].y,cpts[i2].x, cpts[i2].y);
	//exit(42);
	}
    transforms.push_back(t);
    double sumx = orig[i0].x + orig[i1].x + orig[i2].x;
    double sumy = orig[i0].y + orig[i1].y + orig[i2].y;
    centers.push_back(Point(sumx/3.0, sumy/3.0));
    }
return corr;
}

// find the gradients across a set of points in a picture.
void FindGradients(Picture &pb, vector<Point> &pts, double &slopex, double &slopey)
{
double sum_p = 0.0;   // sum of pixel values
double sum_xp = 0.0;  // x times pixel
double sum_x = 0.0;   // sum of the x values
double sum_x2 = 0.0;  // sum of x^2 values
double sum_yp = 0.0;  // y times pixel
double sum_y = 0.0;   // sum of the x values
double sum_y2 = 0.0;  // sum of x^2 values
int N = pts.size();
int w = pb.w;
for(int i=0; i<N; i++) {
    int ix = int(pts[i].x);
    int iy = int(pts[i].y);
    int pixel = pb.DoG[iy*w+ix];
    sum_p += pixel;
    sum_xp += ix*pixel;
    sum_x += ix;
    sum_x2 += ix*ix;
    sum_yp += iy*pixel;
    sum_y += iy;
    sum_y2 += iy*iy;
    }
double avgx = sum_x/N;
slopex = (sum_xp - N*avgx*(sum_p/N))/(sum_x2 - N*avgx*avgx);
double avgy = sum_y/N;
slopey = (sum_yp - N*avgy*(sum_p/N))/(sum_y2 - N*avgy*avgy);
fprintf(of, "Slope bias per 1000 pixels: in x %f, in y %f\n", slopex*1000.0, slopey*1000.0);
}

// write a picture from an normalized array of doubles.
void WriteMonoTiffFromDouble(const char *filename, vector<double> &b, int w, int h)
{
uint8 *buffer = (uint8 *)malloc(b.size()*sizeof(uint8));
for(int i=0; i<b.size(); i++) {
    int pix = int(127 + b[i]*35.0);
    if (pix < 0) pix = 0;
    if (pix > 255) pix = 255;
    buffer[i] = pix;
    }
#ifdef WRITE_DEBUG_IMAGES
WriteMonochromeTIFF(filename, buffer, w, h, 8);
#endif
free(buffer);
}

CD FFT_r2c_lookup(vector<CD> &c, int Nx, int Ny, int lx, int ly)
{
int M = Nx/2+1;
if (lx < 0) lx += Nx;
if (ly < 0) ly += Ny;
if (lx > Nx/2)
    return conj(c[(Nx-lx) + ly*M]);
return c[lx + ly*M];
}
static int NextFFT = 0;
// See if two images match in the Fourier domain
// wvlen is the maximum size of pixels
// Point values have already been normalized
double FourierMatch(vector<Point> &pts, vector<double> &a, vector<double> &b, int wvlen, bool write_images, const char *msg)
{

double pwr = 0.0;
for(int i=0; i<a.size(); i++)
    pwr += (a[i]-b[i])*(a[i]-b[i]);
fprintf(of, "%s: Energy per pixel in diff image %f\n", msg, pwr/a.size());

// step 1 - find the angle that minimizes the area of the bounding rectangle
// Trying all angle with several million points could be slow, so find a bounding set of points.
// We'll take the min and max x of each Y coordinate (called the outline) and transform that
double xmin = BIG, xmax = -BIG, ymin = BIG, ymax = -BIG;
for(int i=0; i<pts.size(); i++) {
    xmin = fmin(xmin, floor(pts[i].x));
    xmax = fmax(xmax, ceil (pts[i].x));
    ymin = fmin(ymin, floor(pts[i].y));
    ymax = fmax(ymax, ceil (pts[i].y));
    }
vector<double> minx(int(ymax-ymin) + 1, xmax + 1);
vector<double> maxx(int(ymax-ymin) + 1, xmin - 1);
for(int i=0; i<pts.size(); i++) {
    int yi = int(floor(pts[i].y-ymin));   // y used to index into the array
    minx[yi] = min(minx[yi], pts[i].x);
    maxx[yi] = max(maxx[yi], pts[i].x);
    }
vector<Point> outline;
for(int i=0; i<pts.size(); i++) {
    int yi = int(floor(pts[i].y-ymin));   // y used to index into the array
    if (minx[yi] < maxx[yi]) {
	outline.push_back(Point(minx[yi], pts[i].y));
	outline.push_back(Point(maxx[yi], pts[i].y));
	}
    }
//
// Finally, now that we've got a smaller set, find the best angle
double best_area = 1.0E30;
double best_angle = 0.0;
for(int angle = -45; angle <= 45; angle++) {
    double rad = angle/180.0*PI;  // angle in radians
    double co = cos(rad), si = sin(rad);
    tform t(co, -si, 0, si, co, 0.0);
    double xmin = BIG, xmax = -BIG, ymin = BIG, ymax = -BIG;
    for(int i=0; i<outline.size(); i++) {
        Point p = outline[i];
        Transform(p,t);
	xmin = fmin(xmin, p.x);
	xmax = fmax(xmax, p.x);
	ymin = fmin(ymin, p.y);
	ymax = fmax(ymax, p.y);
	}
    double area = (xmax-xmin)*(ymax-ymin);
    if (area < best_area) {
	best_area = area;
        best_angle = rad;
        }
    }
fprintf(of, "%s: best angle is %f radians (%f degrees), area %f\n", msg, best_angle, best_angle*180/PI, best_area);

// create the best transform
tform bt(cos(best_angle), -sin(best_angle), 0.0, sin(best_angle), cos(best_angle), 0.0);

xmin = BIG; xmax = -BIG; ymin = BIG; ymax = -BIG;
double xr = -BIG;
double yr = -BIG;
for(int i=0; i<pts.size(); i++) {
    Point p = pts[i];
    Transform(p, bt);
    xmin = fmin(xmin, floor(p.x));
    xmax = fmax(xmax,  ceil(p.x));
    ymin = fmin(ymin, floor(p.y));
    ymax = fmax(ymax,  ceil(p.y));
    xr   = fmax(xr,  p.x);   
    yr   = fmax(yr,  p.y);   
    }
if (xr == xmax) // then the biggest was an exact integer
    xmax += 1;  // not needed, in theory, since then interpolation will add 0 the next (off the array end) element.
                // but in the interest of keeping all accesses within the array, we bump the size in this case
if (yr == ymax) // same for y
    ymax += 1;  
// find the FFT size needed.
int Nx = int(xmax - xmin) + 1;
int Ny = int(ymax - ymin) + 1;
// Make sure they are at least even, for ease in indexing frequencies
if (Nx & 1) Nx++;
if (Ny & 1) Ny++;
fprintf(of, "range x %f %f, y %f %f, so use Nx=%d, Ny=%d\n", xmin, xmax, ymin, ymax, Nx, Ny);

// Make images
vector<double> i1(Nx*Ny, 0.0);
vector<double> i2(Nx*Ny, 0.0);
for(int i=0; i<pts.size(); i++) {
    Point p = pts[i];
    Transform(p, bt);
    int ix = int(floor(p.x - xmin));
    int iy = int(floor(p.y - ymin));
    double alpha = p.x - xmin - ix;
    double beta  = p.y - ymin - iy;
    i1[ix   + Nx*iy    ] += (1-alpha)*(1-beta)*a[i];
    i1[ix+1 + Nx*iy    ] += (  alpha)*(1-beta)*a[i];
    i1[ix   + Nx*(iy+1)] += (1-alpha)*(  beta)*a[i];
    i1[ix+1 + Nx*(iy+1)] += (  alpha)*(  beta)*a[i];
    i2[ix   + Nx*iy    ] += (1-alpha)*(1-beta)*b[i];
    i2[ix+1 + Nx*iy    ] += (  alpha)*(1-beta)*b[i];
    i2[ix   + Nx*(iy+1)] += (1-alpha)*(  beta)*b[i];
    i2[ix+1 + Nx*(iy+1)] += (  alpha)*(  beta)*b[i];
    }

// Create the difference image
vector<double> diff(Nx*Ny);
for(int i=0; i<Nx*Ny; i++)
    diff[i] = i1[i] - i2[i];

if (write_images) {
    char fname[32];
    sprintf(fname,"fft%d-a.tif", NextFFT);
    WriteMonoTiffFromDouble(fname,i1,Nx,Ny);
    sprintf(fname,"fft%d-b.tif", NextFFT);
    WriteMonoTiffFromDouble(fname,i2,Nx,Ny);
    double e1 = 0.0, e2 = 0.0, ed = 0.0;
    for(int i=0; i<Nx*Ny; i++) {
        e1 += i1[i]*i1[i];
        e2 += i2[i]*i2[i];
        ed += diff[i]*diff[i];
        }
    fprintf(of,"energies %f %f %f\n", e1, e2, ed);
    sprintf(fname,"fft%d-d.tif", NextFFT);
    WriteMonoTiffFromDouble(fname,diff,Nx,Ny);
    NextFFT++;
    }
// FFT each one
int M = Ny*(Nx/2+1); // number of complex outputs
vector<CD> i1fft(M);
fftw_plan p;
p = fftw_plan_dft_r2c_2d(Ny, Nx, (double *)(&i1[0]), (double (*)[2])(&i1fft[0]), FFTW_ESTIMATE);
fftw_execute(p); // execute the plan
fftw_destroy_plan(p);
vector<CD> i2fft(M);
p = fftw_plan_dft_r2c_2d(Ny, Nx, (double *)(&i2[0]), (double (*)[2])(&i2fft[0]), FFTW_ESTIMATE);
fftw_execute(p); // execute the plan
vector<CD> dfft(M);
p = fftw_plan_dft_r2c_2d(Ny, Nx, (double *)(&diff[0]), (double (*)[2])(&dfft[0]), FFTW_ESTIMATE);
fftw_execute(p); // execute the plan
fftw_destroy_plan(p);
//
// Compare the lowest few coefficients, corresponding to wavelenths longer than wvlen
double total = 0.0;
double dot = 0.0;
int xlim = Nx/wvlen;
int ylim = Ny/wvlen;
for(int x = -xlim; x <= xlim; x++) {
    for(int y = -ylim; y <= ylim; y++) {
	CD v1 = FFT_r2c_lookup(i1fft, Nx, Ny, x, y);
	CD v2 = FFT_r2c_lookup(i2fft, Nx, Ny, x, y);
        double ang_a = atan2(v1.imag(), v1.real());
        double ang_b = atan2(v2.imag(), v2.real());
        double mag_a = sqrt(v1.real()*v1.real() + v1.imag()*v1.imag());
        double mag_b = sqrt(v2.real()*v2.real() + v2.imag()*v2.imag());
        total += mag_a * mag_b;
        dot += v1.real()*v2.real() + v1.imag()*v2.imag();
	//fprintf(of, "a: %f %f %f %f  b: %f %f %f %f\n", v1.real(), v1.imag(), mag_a, ang_a,
	//v2.real(), v2.imag(), mag_b, ang_b);
	}
    }
fprintf(of,"Metric: %f, energy %f\n", dot/total, total);

// Look at the power spectrum of the differences
// The wavelengths of the harmonics range from 2px to Npx, in both X and Y.  So the biggest
// wavelength contained in the picture is the bigger of nx and ny, in pixels
// The DC terms represent even longer wavelengths.  The shortest wavelength is a wave
// traveling at 45 degrees to the axes, which gives a measured wavelength of 2 px on
// each axis.  Note that 1e30 serves as infinity for our purposes
int MM = Nx/2 + 1; // number of components in X (only slightly more than half since input is real
int maxf = max(Nx, Ny) + 1;  // add one for the zero frequency component
vector<double>pspectrum(maxf, 0.0);
for(int x = 0; x < MM; x++) {
    double wavex = x > 0 ? double(Nx)/x : 1.0e30;
    for(int y = 0; y < Ny; y++) {
        int iy = y > Ny/2 ? Ny-y : y; // frequency reflects about the middle of the array
        double wavey = y > 0 ? double(Ny)/iy : 1.0e30;
        double ratio = wavex/wavey;
        double wave = wavex/sqrt(1.0 + ratio*ratio);
        if (wave > maxf-1)
            wave = maxf-1;
        int iwave = int(floor(wave));
        pspectrum[iwave] += norm(dfft[x + MM*y]);
        }
    }
double psum = 0.0;
for(int i=0; i<maxf; i++)
    psum += pspectrum[i];
double sum = 0.0;
int i;
for(i=0; sum < psum*0.90 && i < maxf; i++) {
    sum += pspectrum[i];
    if (i  > 0 && (i % 10) == 0)
	fprintf(of, "%s: cum. to %d = %f\n", msg, i, sum/psum);
    }
fprintf(of,"%s: Power exceeds half of %f at index %d (subsum %f)\n", msg, psum, i, sum/psum);
return dot/total;
}

double InterpolatePixel(double x, double y, vector<double> &image, int N)
{
int ix = int(x);
int iy = int(y);
double alpha = x - ix;
double beta  = y - iy;
double v = (1-alpha)*(1-beta)*image[ix+N*iy]     + (alpha)*(1-beta)*image[ix+1+N*iy] +
	   (1-alpha)*(  beta)*image[ix+N*(iy+1)] + (alpha)*(  beta)*image[ix+1+N*(iy+1)];
return v;
}

// try another metric - by eye, we say a match is good if it has a lot of yellow points.
// take two normalized images, see what percentage of pixels are yellow.  Try several thresholds for now
double PercentYellow(vector<double> &a, vector<double> &b)
{
int yellow;
for(double thr = 0.05; thr < 0.255; thr += 0.05) {
    int red = 0, green = 0; 
    yellow = 0;
    for(int i=0; i<a.size(); i++) {
        // create pixel values, as we usually do
        int pa = int(a[i]*35 + 127);
        int pb = int(b[i]*35 + 127);
        if (pa < 0) pa = 0; if (pa > 255) pa = 255;
        if (pb < 0) pb = 0; if (pb > 255) pb = 255;
        double ratio = double(pa)/double(pb);
        if (ratio > 1 + thr)
	    red++;
        else if (ratio < 1.0/(1+thr))
	    green++;
	else
	    yellow++;
	}
    fprintf(of, " Thresh %6.1f, red %6.1f   yellow %6.1f   green %6.1f\n", thr*100.0, 
     double(red)/a.size()*100.0, double(yellow)/a.size()*100.0, double(green)/a.size()*100.0);
    }
return double(yellow)/a.size();
}
// Starting with an approximate transform, find detailed correspondence
// Values are returned by setting entries in the array ids (to tell what the regions are)
// and the appending the corresponding entries to map1.transforms.
// The pixel contents in the map 'id' are set to index of transform, plus 10.
// Write triangles to file 'ftri', if not NULL
void RegionToRegionMap(Picture &pa, ConnRegion &a, Picture &pb, ConnRegion &b, 
tform &tr_guess, uint16 *ids, ffmap &map1,
FILE *flog, FILE *ftri)
{

// now we need to find how each of these regions map to the picture below
// if indeed they do at all.
int w = pa.w;		// must be the same for a and b, verified earlier
int h = pa.h;
int sc = pa.scale;
int fullw = pa.w*sc;
int fullh = pa.h*sc;
int npixels = b.pts.size();
fprintf(of," Finding detailed mapping - scale = %d\n", sc);
   {
   vector<CD> cache;
   pa.MakeDoGExist(InnerR, OuterR, cache);
   pb.MakeDoGExist(InnerR, OuterR, cache);
   }
fprintf(of," Done with DoG versions\n");

// Here we create a normalized versions.  Initially, just mean=0 and std=1.  
// Now we subtract out gradients, too.
double slope_x, slope_y;
FindGradients(pb, b.pts, slope_x, slope_y);
MeanStd mm;
for(int i=0; i<npixels; i++) {
    int ix = int(b.pts[i].x);
    int iy = int(b.pts[i].y);
    mm.Element(pb.DoG[iy*w+ix] - ix*slope_x - iy*slope_y);
    }
fprintf(of,"below: %d real pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/(w*h));
double mean2, std2;
mm.Stats(mean2, std2);  // find existing statistics
fprintf(of,"On the target image,  mean= %f and std dev = %f\n", mean2, std2);

// Make a 4Kx4K normalized copy, with image in lower left 2Kx2K, so FFT will work.  
// Background pixels map to zero.
// Also, create a map that just tells whether each point is in the region, or not
vector<double> image2(4096*4096, 0.0);
vector<char> IsIn(w*h,0);
for(int i=0; i<npixels; i++) {
    int x = int(b.pts[i].x);
    int y = int(b.pts[i].y);
    double pix = (pb.DoG[y*w + x] -x*slope_x - y*slope_y - mean2)/std2;
    if (abs(pix) > 2.5)    // outlying pixel value
	pix = 0.0;  // will be set to the mean value (0)
    image2[4096*y+x] = pix;
    IsIn[w*y+x] = 1;
    }

// if needed, make a bigger copy that can be used for full resolution 
//vector<double>image3;
//if (sc > 1) { // then we have a copy with more than 2Kx2K resolution
//    int np3 = npixels*sc*sc;
//    MeanStd m3;
//    for(int i=0; i<np3; i++) {
//	if (below[0].orig[i] > 0)
//	    m3.Element(below[0].orig[i]);
//	}
//    double mean3, std3;
//    m3.Stats(mean3, std3);  // find existing statistics
//    fprintf(of,"Of the full resolution target image,  mean= %f and std dev = %f\n", mean3, std3);
//    image3.push_back(0.0);
//    image3.insert(image3.begin(), 4096*4096-1, 0.0); // is there a better way?
//    int bigw = below[0].w*sc;
//    for(int i=0; i<np3; i++) {
//	int y = i/bigw;
//	int x = i-y*bigw;
//	double pix = below[0].orig[i];
//	if (pix == 0)    // background pixels
//	    pix = mean3;  // will be set to the mean value (0)
//	image3[4096*y+x] = (pix-mean3)/std3;
//	}
//    }

// create a map image.  For every pixel in A, tell the closest control point.  0-9  reserved for 'no data',
// for cases such as folds, off the target image, piece too small, etc. .  
// Otherwise if the pixel has value N, the closest point is N-10
//

fprintf(of,"\n-----Starting detailed region mapping ----\n");
vector<vertex> ControlPoints;
vector<triangle> tris;

// now using the guess, keep only those points that map to the below image.  (this could be off
// by 10 pixels or so, so we'll need to fix it later, but it's a good guess
vector<Point> newpts;
for(int i=0; i<a.pts.size(); i++) {
    Point p(a.pts[i]);
    Transform(p, tr_guess);
    int ix = int(p.x);   // approximate target pixel
    int iy = int(p.y);
    if (ix >= 0 && ix <= w-1 && iy >= 0 && iy <= h-1 && IsIn[iy*w+ix] != 0)
	newpts.push_back(a.pts[i]);
    }
// some quality checks should go here
fprintf(of,"Roughly %d pixels map to the image below\n", newpts.size());
if (newpts.size() < MinMapArea) { // region too small?
    fprintf(of, "STAT: Region too small - %d pixels, limit %.1f.\n", newpts.size(), MinMapArea);
    return;
    }

vector<double> spv;  // create a normalized array of source pixel values
                     // remove any gradients as well
FindGradients(pa, a.pts, slope_x, slope_y); // find gradients using whole region
					    // since the overlapping part might be small, so removing
					    // gradients could be misleading
for(int k=0; k<newpts.size(); k++) {
    int ix = int(newpts[k].x);
    int iy = int(newpts[k].y);
    spv.push_back(pa.DoG[ix + iy*w] - ix*slope_x - iy*slope_y);
    }
Normalize(spv);

vector<vertex> SimplePts;
vector<triangle> SimpleTri;
// Now create a rough ouline of the overlap region
if (CreateOutline(newpts, ControlPoints, tris, SimplePts, SimpleTri)) {
    fprintf(of, "STAT: Triangular mesh failed - %d pixels, limit %.1f.\n", newpts.size(), MinMapArea);
    return;  //overlap region could be too small, for example
    }
vector<tform>transforms;
vector<Point> centers;
double corr = ImproveMesh(newpts, spv, // describes source image
 SimplePts, SimpleTri,        // describes triangles in source
 image2,                                       // describes the target image
 tr_guess, DetailFinalThreshold*0.9,           // approximate transform
 transforms, centers,            // results go here
 flog, "affine");

if (transforms.size() == 0)
    return;
else {
   fprintf(of, "Best affine transform: "); PrintTransform(of, transforms[0]);
   tform inv;
   InvertTrans(inv, transforms[0]);
   fprintf(of, "    Inverse transform: "); PrintTransform(of, inv);
   }
// see if result makes sense in the Fourier domain as well
Match mat;
for(int k=0; k<a.pts.size(); k++) {
    Point pt = a.pts[k];
    Transform(pt, transforms[0] );
    if (pt.x >= 0 && pt.x < w-1 && pt.y >= 0.0 && pt.y < h-1 && IsIn[int(pt.x) + int(pt.y)*w] ) {
        // OK, the point maps.  Fill in the match array
        int ix = int(a.pts[k].x);
        int iy = int(a.pts[k].y);
        mat.pts.push_back(a.pts[k]);
        mat.a.push_back(pa.DoG[ix+iy*w] - ix*slope_x - iy*slope_y);
        mat.b.push_back(InterpolatePixel(pt.x, pt.y, image2, 4096) );
        }
    }
// check it
Normalize(mat.a);
Normalize(mat.b);
{
    double sum = 0.0;
    int N = mat.a.size();
    for(int i=0; i<N; i++)
	sum += mat.a[i] * mat.b[i];
    //for(int harm=1; harm<20; harm++) {
        //double fm = FourierMatch(mat.pts, mat.a, mat.b, harm, harm==1);
        //fprintf(of, "Affine Triangle: %d points, corr %f, %d fm %f \n", N, sum/N, harm, fm);
        //}
    double fm = FourierMatch(mat.pts, mat.a, mat.b, 25, false, "AFF");
    PercentYellow(mat.a, mat.b);
    fprintf(of, "Affine Triangle: %d points, corr %f, Fourier %f\n", N, sum/N, fm );
    if(fm < InitialFourierMetric) {
        fprintf(of  ,"Kicking out - Fourier metric too low (%f), limit %f\n", fm, InitialFourierMetric);
        fprintf(flog,"STAT: Kicking out - Fourier metric too low (%f), limit %f\n", fm, InitialFourierMetric);
        for(int wvlen=5; wvlen <= 40; wvlen += 5) {
            double ffm = FourierMatch(mat.pts, mat.a, mat.b, wvlen, wvlen==5, "AF2");
            fprintf(of, " wavelength %d, metric %f\n", wvlen, ffm);
            }
	return;
	}
    }
//
if (JustOne != 0.0) {
    ControlPoints = SimplePts;
    tris = SimpleTri;
    }
else { // improve the guess here
    tr_guess = transforms[0];
    transforms.clear();
    centers.clear();

    corr = ImproveMesh(newpts, spv, // describes source image
     ControlPoints, tris,        // describes triangles in source
     image2,                                       // describes the target image
     tr_guess, DetailFinalThreshold,               // approximate transform
     transforms, centers,            // results go here
     flog, "deformable mesh");

    if(transforms.size() == 0)
	return;
    }

// Collect values for the final QA step
Match allp;                            // all pixels in the patch
vector<Match> matches(tris.size());    // and indexed by triangle
fprintf(of," Remapping %d points\n", a.pts.size());
for(int k=0; k<a.pts.size(); k++) {
    Point pt = a.pts[k];
    int reg = BestTriangle(pt, tris, ControlPoints);
    Transform(pt, transforms[reg] );
    if (pt.x >= 0 && pt.x < w-1 && pt.y >= 0.0 && pt.y < h-1 && IsIn[int(pt.x) + int(pt.y)*w] ) {
        // OK, the point maps.  Fill in the match array
        int ix = int(a.pts[k].x);
        int iy = int(a.pts[k].y);
        double pva = pa.DoG[ix+iy*w] - ix*slope_x - iy*slope_y;
        double pvb = InterpolatePixel(pt.x, pt.y, image2, 4096);
        matches[reg].pts.push_back(a.pts[k]);
        matches[reg].a.push_back(pva);
        matches[reg].b.push_back(pvb);
        allp.pts.push_back(a.pts[k]);
        allp.a.push_back(pva);
        allp.b.push_back(pvb);
	}
    }

Normalize(allp.a);
Normalize(allp.b);
double dfm = FourierMatch(allp.pts, allp.a, allp.b, 25, false, "DEF");
fprintf(of, "All points, deformable, gives %f\n", dfm);

// check each triangle one by one.  Also accumulate a grand weighted total.
int sum_pts = 0;
double weighted_sum = 0.0;
double weighted_yellow = 0.0;
for(int k=0; k<tris.size(); k++) {
    Normalize(matches[k].a);
    Normalize(matches[k].b);
    double sum = 0.0;
    int N = matches[k].a.size();
    for(int i=0; i<N; i++)
	sum += matches[k].a[i] * matches[k].b[i];
    double fm = FourierMatch(matches[k].pts, matches[k].a, matches[k].b, 25, false, "TRI");
    double y = PercentYellow(matches[k].a, matches[k].b);
    fprintf(of, "Triangle %d, %d points, corr %f, fm %f \n", k, N, sum/N, fm);
    sum_pts += N;
    weighted_sum += fm*N;
    weighted_yellow += y*N;
    }
double score = weighted_sum/sum_pts;
double yell = weighted_yellow/sum_pts;
fprintf(of,"STAT: Overall %d points, %f corr, dfm %f, weighed Fourier metric %f, weighted yellow %6.4f\n", 
 sum_pts, corr, dfm, score, yell);
fprintf(flog,"STAT: Overall %d points, %f corr, dfm %f, weighed Fourier metric %f, weighted yellow %6.4f\n", 
 sum_pts, corr, dfm, score, yell);
if (dfm < FinalFourierMetric || yell < FinalYellowLimit) {
    fprintf(of, "STAT: rejected\n");
    return;
    }

// OK, we are finally convinced the data is OK
// copy the transforms and centers to the proper place
int next_id = map1.transforms.size() + 10;
for(int k=0; k<tris.size(); k++) {
    map1.transforms.push_back(transforms[k]);
    map1.centers.push_back(centers[k]);
    }

// Now, using the final transforms, figure out which points map to the below image, and fill in the map
fprintf(of," Final Remapping %d points\n", a.pts.size());
for(int k=0; k<a.pts.size(); k++) {
    Point pt = a.pts[k];
    int reg = BestTriangle(pt, tris, ControlPoints);
    Transform(pt, map1.transforms.at(next_id-10+reg) );
    if (pt.x >= 0 && pt.x < w-1 && pt.y >= 0.0 && pt.y < h-1 && IsIn[int(pt.x) + int(pt.y)*w] ) {
        // OK, the point maps.  Fill in the match array
        int ix = int(a.pts[k].x);
        int iy = int(a.pts[k].y);
	pt = a.pts[k]; // get the original point again
	for(int x=0; x<sc; x++)  // array ids is at original scale
	    for(int y=0; y<sc; y++)
		ids[sc*int(pt.x)+x + fullw*(sc*int(pt.y)+y)] = next_id + reg;
	}
    }

// One final print, to help in debugging
for(int i=next_id-10; i<map1.centers.size(); i++) {
    Point pt = map1.centers[i];
    fprintf(of," center %f %f :", pt.x, pt.y);
    PrintTransform(of, map1.transforms[i]);
    Point p2 = pt;
    Transform(p2, map1.transforms[i]);
    fprintf(of, "Mapping region %d xy= %f %f to region %d xy= %f %f\n", a.id, pt.x, pt.y, b.id, p2.x, p2.y);
    }

// And print out a mapping of the triangles, for trakEM
if (ftri != NULL) {
    for(int j=0; j<tris.size(); j++) {
        Point ps[3];
        for(int k=0; k<3; k++)
            ps[k] = Point(ControlPoints[tris[j].v[k]].x, ControlPoints[tris[j].v[k]].y);
        fprintf(ftri, "%f %f %f %f %f %f       ", ps[0].x, ps[0].y, ps[1].x, ps[1].y, ps[2].x, ps[2].y);
        for(int k=0; k<3; k++)
            Transform(ps[k], transforms[j]);
        fprintf(ftri, "%f %f %f %f %f %f\n", ps[0].x, ps[0].y, ps[1].x, ps[1].y, ps[2].x, ps[2].y);
        }
    }
}

// write out the parameters we use during analysis
void PrintParams(FILE *of) 
{
// First print out the values of the parameters used:
fprintf(of, "DIT  = %10.3f; When we call the detailed optimizer, starting correlation should be this good at least\n", DetailInitialThreshold);
fprintf(of, "DFT  = %10.3f; After we call the detailed optimizer, correlation should be this good at least\n", DetailFinalThreshold);
fprintf(of, "TAC  = %10.3f; Two transforms are 'close' if within this value (at the 2K scale)\n", TransAreClose);
fprintf(of, "MNL  = %10.3f; Triangle edges should be between minl and 2*minl;\n", minl);
fprintf(of, "MTA  = %10.3f; Minimum allowable size of a triangle (pixels)\n", MinTriangleArea);
fprintf(of, "MNC  = %10.3f; Minimum normalized correlation to define a crude match\n", MinNormCorrelation);
fprintf(of, "MNG  = %10.3f; Minimum normalized correlation to define a great match (can stop here)\n", MinNormGreat);
fprintf(of, "MMA  = %10.3f; This much area should map from one region to another\n", MinMapArea);
fprintf(of, "NCP  = %10.3f; Normalized corr not computed on sizes smaller than this\n", NumCorrPatchSize);
fprintf(of, "SCALE= %10.3f; Scale of below image used in initial search\n", ApproxScale);
fprintf(of, "XSCALE= %10.3f; X only scale of below image used in initial search\n", InitXScale);
fprintf(of, "YSCALE= %10.3f; Y only scale of below image used in initial search\n", InitYScale);
fprintf(of, "SKEW = %10.3f; Skew the image by this much before matching\n", InitialSkew);
fprintf(of, "ONE  = %10.3f; If non-zero, use just one affine transform\n", JustOne);
fprintf(of, "INR  = %10.3f; Inner radius of difference of gaussians, in pixels\n", InnerR);
fprintf(of, "OUTR = %10.3f; Outer radius of difference of gaussians, in pixels (0 -> not used)\n", OuterR);
fprintf(of, "IFM =  %10.3f; Initial Fourier Metric must be this good after affine match\n", InitialFourierMetric);
fprintf(of, "FFM  = %10.3f; Final Fourier Metric must be this good after deformable mesh\n", FinalFourierMetric);
fprintf(of, "CTR  = %10.3f; Center angle for search\n", AngleCenter);
fprintf(of, "FYL  = %10.3f; This fraction of pixels must be yellow\n", FinalYellowLimit);
}


// The main routine the pipeline should call.
void PipelineDeformableMap(int w, int h, // size of all images
 uint8 *AboveRaster,                // the higher layer
 uint8 *fold_mask_above,            // 0 for fold, 1=first connected region, 2 = second, etc.
 uint8 *BelowRaster,                // the lower layer
 uint8 *fold_mask_below,
 uint16 *map_mask,                  // the resulting map.  <10 means no mapping, 10+i means use transform i
 int &Ntrans,                       // how many tranforms result?  (returned value)
 double * &array_of_transforms,     // array of these values.
 FILE *fout,			    // Most of the output goes here
 FILE *flog)                        // A few detailed things get written here
{
of = fout;
PrintParams(fout);
PrintParams(flog);

Picture above, below; // convert to this for sake of compatibility
above.original = AboveRaster;
below.original = BelowRaster;
above.w = below.w = w;
above.h = below.h = h;

//create smaller (at most 2Kx2K pixels, if original was bigger
above.DownsampleIfNeeded();
below.DownsampleIfNeeded();

// Create the connected region lists.  Note that the connected region lists are always at the reduced resolution, if this is used.
vector<ConnRegion> Acr, Bcr;
FindConnRegions(Acr, fold_mask_above, w, h, above.scale);
FindConnRegions(Bcr, fold_mask_below, w, h, below.scale);

for(int j=0; j<w*h; j++)
    map_mask[j] = 0;
ffmap map1;  // data structure holding feature to feature mapping
FILE *ftri = fopen("Triangles.txt","w");
// now find all the mappings from each connected region to each connected region
for(int i = 0; i <Acr.size(); i++) {
    for(int j=0; j < Bcr.size(); j++) {
        fprintf(of,"Looking for mapping from region %d of 'above' to region '%d' of below\n", i, j);
        vector<tform> guesses;
        if (ApproximateMatch(above, Acr[i], below, Bcr[j], guesses, flog)) {
            int num = 0; // count how many approximate transforms result in good matches
	    for(int k=0; k<guesses.size(); k++) {
                fprintf(of,"STAT: try %d\n", k);
                int before = map1.transforms.size();
                RegionToRegionMap(above, Acr[i], below, Bcr[j], guesses[k], map_mask, map1, flog, ftri);
		num += (map1.transforms.size() != before);
                if (map1.transforms.size() != before)
		    break;  // BLOOT - this makes next warning useless
		}
            if (num > 1)
		fprintf(of, "Warning!  %d approximate transforms matched!\n", num);
	    }
	}
    }
if (ftri != NULL)
    fclose(ftri);
// Plausibility check
Ntrans = map1.transforms.size();
if (Ntrans > 1) {
    double ang_min = PI;
    double ang_max = -PI;
    for(int i=0; i<Ntrans; i++) {
	double ang = atan2(map1.transforms[i].t[1], map1.transforms[i].t[0]);
        ang_min = fmin(ang_min, ang);
        ang_max = fmax(ang_max, ang);
	}
    fprintf(fout, "Angle span: min, max, delta = %f %f %f\n", ang_min, ang_max, ang_max-ang_min);
    if (ang_max -ang_min > 0.2) {
	fprintf(fout, "STAT: Angle span TOO BIG: %f %f %f\n", ang_min, ang_max, ang_max-ang_min);
	fprintf(flog, "STAT: Angle span TOO BIG: %f %f %f\n", ang_min, ang_max, ang_max-ang_min);
	}
    }

// Now package the transforms in the way the MatLab code would like to see them.
fprintf(of,"Modifying transforms: above.scale=%d\n", above.scale);
fprintf(of  ,"Returning %d transforms\n", Ntrans);
fprintf(flog,"Returning %d transforms\n", Ntrans);
array_of_transforms = (double *) malloc(Ntrans*6*sizeof(double));
for(int i=0; i<Ntrans; i++) {
    map1.transforms[i].t[2] *= above.scale;
    map1.transforms[i].t[5] *= above.scale;
    map1.transforms[i].ToMatlab(); // re-order the components
    for(int j=0; j<6; j++) {
        array_of_transforms[i*6+j] = map1.transforms[i].t[j];  // maybe need to re-order here
	}
    }

// write out DoG images for debugging
int npixels = above.w*above.h;

// If the DoG images exist, write them out for debugging
if (above.DoG.size() > 0) {
    //WriteMonochromeTIFF("DoGa.tif", &above.DoG[0], w, h, 8); // 8 bits per pixel
    }
if (below.DoG.size() > 0) {
    //WriteMonochromeTIFF("DoGb.tif", &below.DoG[0], w, h, 8); //  bits per pixel
    }
}
bool operator==(const vertex &a, const vertex &b){return a.x == b.x && a.y == b.y;}
bool operator< (const vertex &a, const vertex &b){return a.x < b.x || (a.x == b.x && a.y < b.y);}

void MakeLambda(vector<double> &l, Point p, vector<Point> &cpts, int i0, int i1, int i2, int i3)
{
// 1 3
// 0 2
double alpha = (p.x-cpts[i0].x)/(cpts[i2].x-cpts[i0].x);
double  beta = (p.y-cpts[i0].y)/(cpts[i1].y-cpts[i0].y);
l[i0] = (1-alpha)*(1-beta);
l[i1] = (1-alpha)*beta;
l[i2] =    alpha*(1-beta);
l[i3] =    alpha*beta;
}

// Tells if a point is inside a given rectangle
bool Inside(int x, int y, int x1, int y1, int x2, int y2)
{return x > x1 && x < x2 && y > y1 && y < y2;}

Point PointAvg(vector<Point> &pts, int i1, int i2, int i3, int i4)
{
Point P( (pts[i1].x+pts[i2].x+pts[i3].x+pts[i4].x)/4, (pts[i1].y+pts[i2].y+pts[i3].y+pts[i4].y)/4 );
return P;
}

// This is the function that defines a legal overlap for this application.  nx and ny are the sizes of the overlap region.
bool InSectionLegal(int nx, int ny, void *a)
{
double n = nx * ny;
return n > 30000 && (nx >= 1500 || ny > 1500);  // 30,000 total pixels, and at least one side is 1500 pixels long
}

//Here we have to code for aligning two images on the same section.  The basic idea is that they are are overlapping
//at the edge by some unknown amount.  Find a few exact correspondence points.
//
// Tries to find a mapping from one connected region to another.
// a = above, b = below, we look for a mapping from a->b
// The picture has the whole image, the connected region tells what part we are concerned with
void InSectionOverlap(int w, int h, uint8 *pa, uint8 *pb, int &Npts, double * &apts, double * &bpts, FILE *flog)
{
of = flog;
Picture a, b; // convert to this for sake of compatibility
a.original = pa;
b.original = pb;
a.w = b.w = w;
a.h = b.h = h;

//create smaller (at most 2Kx2K pixels, if original was bigger
a.DownsampleIfNeeded();
b.DownsampleIfNeeded();
w = a.w;
h = a.h;
// now we need to find how each of these regions overlap at the edges, if indeed they do at all.

// First, for image b, compute an image with just the frame
const double frame = 0.10 ;  // size of the frame
int x1 = int(frame*w);  // coordinates of hole in the center
int x2 = w-1-x1;
int y1 = int(frame*h);
int y2 = h-1-y1;
fprintf(of, "w %d, h %d, x1 %d, y1 %d, x2, %d, y2 %d\n", w, h, x1, y1, x2,y2);
MeanStd mm;
for(int ix=0; ix<w; ix++) {
    for(int iy=0; iy<h; iy++) {
	if (!Inside(ix,iy,x1,y1,x2,y2))
            mm.Element(b.raster[iy*w+ix]);
        }
    }
fprintf(of,"below frame: %d pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/(w*h));
double meanb, stdb;
mm.Stats(meanb, stdb);  // find existing statistics
fprintf(of,"On the target frame,  mean= %f and std dev = %f\n", meanb, stdb);

vector<double>image2(4096*4096, 0.0);
vector<Point> bp;
vector<double> bv;
for(int ix=0; ix<w; ix++) {
    for(int iy=0; iy<h; iy++) {
	if (!Inside(ix,iy,x1,y1,x2,y2)) {
            bp.push_back(Point(ix,iy));
            double pix = (b.raster[iy*w + ix] - meanb)/stdb;
            bv.push_back(pix);
            image2[iy*4096+ix] = pix;
	    }
	}
    }

// Now, for picture A, create a point list, with source pixel values.
mm.Reset();
for(int ix=0; ix<w; ix++) {
    for(int iy=0; iy<h; iy++) {
	if (!Inside(ix,iy,x1,y1,x2,y2))
            mm.Element(a.raster[iy*w+ix]);
        }
    }
fprintf(of,"above frame: %d pixels, %f percent\n", mm.HowMany(), mm.HowMany()*100.0/(w*h));
double meana, stda;
mm.Stats(meana, stda);  // find existing statistics
fprintf(of,"On the source frame,  mean= %f and std dev = %f\n", meana, stda);

vector<Point> pts;
vector<double> vals;
for(int ix=0; ix<w; ix++) {
    for(int iy=0; iy<h; iy++) {
	if (!Inside(ix,iy,x1,y1,x2,y2)) {
	    pts.push_back(Point(ix,iy));
            vals.push_back(a.raster[iy*w+ix]);
	    }
	}
    }
Normalize(vals);

// Now call the FFT routine to do the actual work....
double dx, dy;
vector<CD> ftc;
double c = FindNormCorrelation(pts, vals, bp, bv, dx, dy, 0, 0, 4000, of, InSectionLegal, NULL, NULL, NULL, ftc);
if (c < MinNormCorrelation) {  // empirical
    fprintf(of, "Can't find a good starting point - maximum c = %f\n", c);
    Npts = 0;
    return;
    }
fprintf(of, "dx, dy = %f %f\n", dx, dy);
int idx = int(floor(dx+0.5));   // round not compatible between LINUX and WINDOWS, apparently
int idy = int(floor(dy+0.5));
int xl1, yb1, xl2, yb2;
int xr1, yt1, xr2, yt2;
BoxesFromShifts(w, w, w, w, idx, idy, xl1, yb1, xr1, yt1, xl2, yb2, xr2, yt2);

// here we assume the box is horizontal.
vector<Point> cpts(8,Point(0.0,0.0));
vector<vector<double> >lambdas;
vector<double>spv;  // source pixel values
bool horizontal = xr1-xl1 > yt1-yb1;
double xm1 = xl1 + 0.25*(xr1-xl1);
double xm2 = xl1 + 0.75*(xr1-xl1);
double ym1 = yb1 + 0.25*(yt1-yb1);
double ym2 = yb1 + 0.75*(yt1-yb1);
if (horizontal) {
    // Arrange the points like this:
    // 1    3         5    7
    // 0    2         4    6
    cpts[1] = Point(xl1, yt1); cpts[3] = Point(xm1, yt1); cpts[5] = Point(xm2, yt1); cpts[7] = Point(xr1,yt1);
    cpts[0] = Point(xl1, yb1); cpts[2] = Point(xm1, yb1); cpts[4] = Point(xm2, yb1); cpts[6] = Point(xr1,yb1);
    }
else {  // box is more vertical
    cpts[6] = Point(xl1, yt1); cpts[7] = Point(xr1, yt1);    // 6 7
    cpts[4] = Point(xl1, ym2); cpts[5] = Point(xr1, ym2);    // 4 5
    cpts[2] = Point(xl1, ym1); cpts[3] = Point(xr1, ym1);    // 2 3
    cpts[0] = Point(xl1, yb1); cpts[1] = Point(xr1, yb1);    // 0 1
    }
for(int iy = yb1; iy <= yt1; iy++) {
    for(int ix = xl1; ix <= xr1; ix++) {
        spv.push_back(a.raster[ix + iy*w]);
        vector<double> l(8, 0.0);
        Point p(ix,iy);
        if (horizontal) {
	    if (ix < xm1)
		 MakeLambda(l, p, cpts, 0, 1, 2, 3);
	    else if (ix < xm2)
		 MakeLambda(l, p,  cpts, 2, 3, 4, 5);
	    else
		 MakeLambda(l, p, cpts, 4, 5, 6, 7);
            }
        else {  // box is vertical
	    if (iy < ym1)
		 MakeLambda(l, p, cpts, 0, 2, 1, 3);
	    else if (iy < ym2)
		 MakeLambda(l, p,  cpts, 2, 4, 3, 5);
	    else
		 MakeLambda(l, p, cpts, 4, 6, 5, 7);
            }
        lambdas.push_back(l);
	}
    }
vector<Point> before = cpts;  // start with 8 control points
before.push_back(PointAvg(cpts, 0, 1, 2, 3));  // and add the rectangle centers
before.push_back(PointAvg(cpts, 2, 3, 4, 5));
before.push_back(PointAvg(cpts, 4, 5, 6, 7));
// Now transfer the control points to image b's coordinate system.
for (int i=0; i<8; i++) {
    cpts[i].x += dx;
    cpts[i].y += dy;
    }
// Now improve the correlation by gradient descent
double threshold = 0.5;
c = ImproveControlPts(cpts, lambdas, spv, image2, 4096, flog, "in-section", threshold);
if (c < threshold) {
    fprintf(of, "Correlation not good enough - %f\n", c);
    Npts = 0;
    return;
    }
vector<Point> after = cpts;
after.push_back(PointAvg(cpts, 0, 1, 2, 3));
after.push_back(PointAvg(cpts, 2, 3, 4, 5));
after.push_back(PointAvg(cpts, 4, 5, 6, 7));
// Keep only the points that map into both images.  Source points will always be inside, by construction, but
// some of the target points may fall outside after detailed matching (rotation, for example).
for (int i=0; i<after.size();) {
    if (after[i].x >= 0 && after[i].x <= w-1 && after[i].y >= 0 && after[i].y <= w-1) // it's OK
	i++;
    else {  // image 
	before.erase(before.begin()+i);
	after. erase( after.begin()+i);
        }
    }
Npts = after.size();   
fprintf(of, " %d matching points\n", Npts);
apts = (double *)malloc(2*Npts*sizeof(double));
bpts = (double *)malloc(2*Npts*sizeof(double));

// Need to scale back to original image size, if bigger than 2K
int sc = a.scale;
for(int i=0; i<Npts; i++) {
    fprintf(of, "Point (%7.2f %7.2f) in image 1 maps to point (%7.2f %7.2f) in image 2, delta (%7.2f %7.2f)\n", 
     sc*before[i].x, sc*before[i].y, sc*after[i].x, sc*after[i].y,
     sc*before[i].x-sc*after[i].x,   sc*before[i].y-sc*after[i].y);
    apts[2*i] = sc*before[i].x; apts[2*i+1] = sc*before[i].y;
    bpts[2*i] = sc*after[i].x; bpts[2*i+1] = sc*after[i].y;
    }
return;
}
