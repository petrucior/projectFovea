/* Original code has been submitted by Liu Liu. Here is the copyright.
----------------------------------------------------------------------------------
 * An OpenCV Implementation of SURF
 * Further Information Refer to "SURF: Speed-Up Robust Feature"
 * Author: Liu Liu
 * liuliu.1987+opencv@gmail.com
 *
 * There are still serveral lacks for this experimental implementation:
 * 1.The interpolation of sub-pixel mentioned in article was not implemented yet;
 * 2.A comparision with original libSurf.so shows that the hessian detector is not a 100% match to their implementation;
 * 3.Due to above reasons, I recommanded the original one for study and reuse;
 *
 * However, the speed of this implementation is something comparable to original one.
 *
 * Copyright© 2008, Liu Liu All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or
 * without modification, are permitted provided that the following
 * conditions are met:
 *  Redistributions of source code must retain the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer.
 *  Redistributions in binary form must reproduce the above
 *  copyright notice, this list of conditions and the following
 *  disclaimer in the documentation and/or other materials
 *  provided with the distribution.
 *  The name of Contributor may not be used to endorse or
 *  promote products derived from this software without
 *  specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 */

/*
   The following changes have been made, comparing to the original contribution:
   1. A lot of small optimizations, less memory allocations, got rid of global buffers
   2. Reversed order of cvGetQuadrangleSubPix and cvResize calls; probably less accurate, but much faster
   3. The descriptor computing part (which is most expensive) is threaded using OpenMP
   (subpixel-accurate keypoint localization and scale estimation are still TBD)
*/

/*
KeyPoint position and scale interpolation has been implemented as described in
the Brown and Lowe paper cited by the SURF paper.

The sampling step along the x and y axes of the image for the determinant of the
Hessian is now the same for each layer in an octave. While this increases the
computation time, it ensures that a true 3x3x3 neighbourhood exists, with
samples calculated at the same position in the layers above and below. This
results in improved maxima detection and non-maxima suppression, and I think it
is consistent with the description in the SURF paper.

The wavelet size sampling interval has also been made consistent. The wavelet
size at the first layer of the first octave is now 9 instead of 7. Along with
regular position sampling steps, this makes location and scale interpolation
easy. I think this is consistent with the SURF paper and original
implementation.

The scaling of the wavelet parameters has been fixed to ensure that the patterns
are symmetric around the centre. Previously the truncation caused by integer
division in the scaling ratio caused a bias towards the top left of the wavelet,
resulting in inconsistent keypoint positions.

The matrices for the determinant and trace of the Hessian are now reused in each
octave.

The extraction of the patch of pixels surrounding a keypoint used to build a
descriptor has been simplified.

KeyPoint descriptor normalisation has been changed from normalising each 4x4
cell (resulting in a descriptor of magnitude 16) to normalising the entire
descriptor to magnitude 1.

The default number of octaves has been increased from 3 to 4 to match the
original SURF binary default. The increase in computation time is minimal since
the higher octaves are sampled sparsely.

The default number of layers per octave has been reduced from 3 to 2, to prevent
redundant calculation of similar sizes in consecutive octaves.  This decreases
computation time. The number of features extracted may be less, however the
additional features were mostly redundant.

The radius of the circle of gradient samples used to assign an orientation has
been increased from 4 to 6 to match the description in the SURF paper. This is
now defined by ORI_RADIUS, and could be made into a parameter.

The size of the sliding window used in orientation assignment has been reduced
from 120 to 60 degrees to match the description in the SURF paper. This is now
defined by ORI_WIN, and could be made into a parameter.

Other options like  HAAR_SIZE0, HAAR_SIZE_INC, SAMPLE_STEP0, ORI_SEARCH_INC,
ORI_SIGMA and DESC_SIGMA have been separated from the code and documented.
These could also be made into parameters.

Modifications by Ian Mahon

*/

#ifndef FOVEATED_HESSIANDETECTOR
#define FOVEATED_HESSIANDETECTOR

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"

using namespace std;
using namespace cv;

static const int   SURF_ORI_SEARCH_INC = 5;
static const float SURF_ORI_SIGMA      = 2.5f;
static const float SURF_DESC_SIGMA     = 3.3f;

// Wavelet size at first layer of first octave.
static const int SURF_HAAR_SIZE0 = 9;

// Wavelet size increment between layers. This should be an even number,
// such that the wavelet sizes in an octave are either all even or all odd.
// This ensures that when looking for the neighbours of a sample, the layers
// above and below are aligned correctly.
static const int SURF_HAAR_SIZE_INC = 6;

struct SurfHF
{
    int p0, p1, p2, p3;
    float w;

    SurfHF(): p0(0), p1(0), p2(0), p3(0), w(0) {}
};

inline float calcHaarPattern( const int* origin, const SurfHF* f, int n )
{
    double d = 0;
    for( int k = 0; k < n; k++ )
        d += (origin[f[k].p0] + origin[f[k].p3] - origin[f[k].p1] - origin[f[k].p2])*f[k].w;
    return (float)d;
}

static void
resizeHaarPattern( const int src[][5], SurfHF* dst, int n, int oldSize, int newSize, int widthStep )
{
    float ratio = (float)newSize/oldSize;
    for( int k = 0; k < n; k++ )
    {
        int dx1 = cvRound( ratio*src[k][0] );
        int dy1 = cvRound( ratio*src[k][1] );
        int dx2 = cvRound( ratio*src[k][2] );
        int dy2 = cvRound( ratio*src[k][3] );
        dst[k].p0 = dy1*widthStep + dx1;
        dst[k].p1 = dy2*widthStep + dx1;
        dst[k].p2 = dy1*widthStep + dx2;
        dst[k].p3 = dy2*widthStep + dx2;
        dst[k].w = src[k][4]/((float)(dx2-dx1)*(dy2-dy1));
    }
}

/*
 * Calculate the determinant and trace of the Hessian for a layer of the
 * scale-space pyramid
 */
template <typename T, typename K>
static void calcLayerDetAndTrace( const Mat& sum, int size, int sampleStep,
				  Mat& det, Mat& trace, K fovea, int marginH, int foveaLevel ){
  cout << "calc" << endl;
  const int NX=3, NY=3, NXY=4;
  const int dx_s[NX][5] = { {0, 2, 3, 7, 1}, {3, 2, 6, 7, -2}, {6, 2, 9, 7, 1} };
  const int dy_s[NY][5] = { {2, 0, 7, 3, 1}, {2, 3, 7, 6, -2}, {2, 6, 7, 9, 1} };
  const int dxy_s[NXY][5] = { {1, 1, 4, 4, 1}, {5, 1, 8, 4, -1}, {1, 5, 4, 8, -1}, {5, 5, 8, 8, 1} };

  //foveated parameters
  int k = foveaLevel;

  vector< T > modelParameters = fovea.getParameters();
  Level< T > levelFovea = fovea.getLevelFromFovea( k );
  vector< T > params = levelFovea.boundingBox();

  int deltax = params[0].x;
  int deltay = params[0].y;
  int skx = params[1].x;
  int sky = params[1].y;

  //margin ref: centro da wavelet
  //margin_x ref: centro da wavelet
  int margin_x = MAX(marginH, deltax);
  int margin_y = MAX(marginH, deltay);

  //limit_x ref: centro da wavelet
  int limit_x = MIN(deltax + skx, modelParameters[2].x - marginH);
  int limit_y = MIN(deltay + sky, modelParameters[2].y - marginH);

  //sum_i ref: comeco da wavelet
  int sum_i, sum_j;
  sum_i = margin_y - size/2;

  //DEBUG
  /*
    std::cout << "Computando a imagem Hessiana" << std::endl;
    std::cout << "Margin H = " << marginH << std::endl;
    std::cout << "foveaLevel = " << foveaLevel << std::endl;
    std::cout << "fovea = " << modelParameters[3].x << " " << modelParameters[3].y << std::endl;
    std::cout << "delta = " << deltax << " " << deltay << std::endl;
    std::cout << "A wavelet vai de " << margin_x << " até " << limit_x << std::endl;
    std::cout << "Pulando de " << sampleStep << " em " << sampleStep << std::endl;
  */
  
  SurfHF Dx[NX], Dy[NY], Dxy[NXY];
  
  if( size > sum.rows-1 || size > sum.cols-1 )
    return;
  
  resizeHaarPattern( dx_s , Dx , NX , 9, size, sum.cols );
  resizeHaarPattern( dy_s , Dy , NY , 9, size, sum.cols );
  resizeHaarPattern( dxy_s, Dxy, NXY, 9, size, sum.cols );

  for(int i = 0; sum_i + size/2 <= limit_y; i++, sum_i += sampleStep ) {
    sum_j = margin_x - size/2;
    const int* sum_ptr = sum.ptr<int>(sum_i, sum_j);
    float* det_ptr = &det.at<float>(i, 0);
    float* trace_ptr = &trace.at<float>(i, 0);
    for(int j = 0; sum_j + size/2 <= limit_x; sum_j += sampleStep, j++ ) {
      float dx  = calcHaarPattern( sum_ptr, Dx , 3 );
      float dy  = calcHaarPattern( sum_ptr, Dy , 3 );
      float dxy = calcHaarPattern( sum_ptr, Dxy, 4 );
      sum_ptr += sampleStep;
      det_ptr[j] = dx*dy - 0.81f*dxy*dxy;
      trace_ptr[j] = dx + dy;
    }
  }
}


/*
 * Maxima location interpolation as described in "Invariant Features from
 * Interest Point Groups" by Matthew Brown and David Lowe. This is performed by
 * fitting a 3D quadratic to a set of neighbouring samples.
 *
 * The gradient vector and Hessian matrix at the initial keypoint location are
 * approximated using central differences. The linear system Ax = b is then
 * solved, where A is the Hessian, b is the negative gradient, and x is the
 * offset of the interpolated maxima coordinates from the initial estimate.
 * This is equivalent to an iteration of Netwon's optimisation algorithm.
 *
 * N9 contains the samples in the 3x3x3 neighbourhood of the maxima
 * dx is the sampling step in x
 * dy is the sampling step in y
 * ds is the sampling step in size
 * point contains the keypoint coordinates and scale to be modified
 *
 * Return value is 1 if interpolation was successful, 0 on failure.
 */
static int
interpolateKeypoint( float N9[3][9], int dx, int dy, int ds, KeyPoint& kpt )
{
    Vec3f b(-(N9[1][5]-N9[1][3])/2,  // Negative 1st deriv with respect to x
            -(N9[1][7]-N9[1][1])/2,  // Negative 1st deriv with respect to y
            -(N9[2][4]-N9[0][4])/2); // Negative 1st deriv with respect to s

    Matx33f A(
        N9[1][3]-2*N9[1][4]+N9[1][5],            // 2nd deriv x, x
        (N9[1][8]-N9[1][6]-N9[1][2]+N9[1][0])/4, // 2nd deriv x, y
        (N9[2][5]-N9[2][3]-N9[0][5]+N9[0][3])/4, // 2nd deriv x, s
        (N9[1][8]-N9[1][6]-N9[1][2]+N9[1][0])/4, // 2nd deriv x, y
        N9[1][1]-2*N9[1][4]+N9[1][7],            // 2nd deriv y, y
        (N9[2][7]-N9[2][1]-N9[0][7]+N9[0][1])/4, // 2nd deriv y, s
        (N9[2][5]-N9[2][3]-N9[0][5]+N9[0][3])/4, // 2nd deriv x, s
        (N9[2][7]-N9[2][1]-N9[0][7]+N9[0][1])/4, // 2nd deriv y, s
        N9[0][4]-2*N9[1][4]+N9[2][4]);           // 2nd deriv s, s

    Vec3f x = A.solve(b, DECOMP_LU);

    bool ok = (x[0] != 0 || x[1] != 0 || x[2] != 0) &&
        std::abs(x[0]) <= 1 && std::abs(x[1]) <= 1 && std::abs(x[2]) <= 1;

    if( ok )
    {
        kpt.pt.x += x[0]*dx;
        kpt.pt.y += x[1]*dy;
        kpt.size = (float)cvRound( kpt.size + x[2]*ds );
    }
    return ok;
}

// Multi-threaded construction of the scale-space pyramid
template <typename T, typename K>
struct SURFBuildInvoker : ParallelLoopBody {

  SURFBuildInvoker( const Mat& _sum, const std::vector<int>& _sizes,
		    const std::vector<int>& _sampleSteps,
		    std::vector<Mat>& _dets, std::vector<Mat>& _traces,
		    K _fovea, vector<int>& _margin, vector<int>& _foveaLevel){
    sum = &_sum;
    sizes = &_sizes;
    sampleSteps = &_sampleSteps;
    dets = &_dets;
    traces = &_traces;
    fovea = &_fovea;
    margin = &_margin;
    foveaLevel = &_foveaLevel;
  }
  
  void operator()(const Range& range) const CV_OVERRIDE {
    for( int i=range.start; i<range.end; i++ ){
      if((*foveaLevel)[i] == -1) continue;
      calcLayerDetAndTrace< T, K >( *sum, (*sizes)[i], (*sampleSteps)[i], (*dets)[i], (*traces)[i], *fovea, (*margin)[i], (*foveaLevel)[i] );
    }
  }
  
  const Mat *sum;
  const std::vector<int> *sizes;
  const std::vector<int> *sampleSteps;
  const vector<int> *foveaLevel;
  const vector<int> *margin;
  std::vector<Mat>* dets;
  std::vector<Mat>* traces;
  K* fovea;
  
};

// Multi-threaded search of the scale-space pyramid for keypoints
template <typename T, typename K>
struct SURFFindInvoker : ParallelLoopBody {
  SURFFindInvoker( const Mat& _sum, const Mat& _mask_sum,
		   const std::vector<Mat>& _dets, const std::vector<Mat>& _traces,
		   const std::vector<int>& _sizes, const std::vector<int>& _sampleSteps,
		   const std::vector<int>& _middleIndices, std::vector<KeyPoint>& _keypoints,
		   int _nOctaveLayers, float _hessianThreshold, K _fovea,
		   vector<int>& _margin, vector<int>& _foveaLevel ) {
    sum = &_sum;
    mask_sum = &_mask_sum;
    dets = &_dets;
    traces = &_traces;
    sizes = &_sizes;
    sampleSteps = &_sampleSteps;
    middleIndices = &_middleIndices;
    keypoints = &_keypoints;
    nOctaveLayers = _nOctaveLayers;
    hessianThreshold = _hessianThreshold;
    fovea = &_fovea;
    margin = &_margin;
    foveaLevel = &_foveaLevel;
  }
  
  static void findMaximaInLayer( const Mat& sum, const Mat& mask_sum,
				 const std::vector<Mat>& dets, const std::vector<Mat>& traces,
				 const std::vector<int>& sizes, std::vector<KeyPoint>& keypoints,
				 int octave, int layer, float hessianThreshold, int sampleStep,
				 K fovea, int marginH, int foveaLevel );
  
  void operator()(const Range& range) const {
    for( int i=range.start; i<range.end; i++ ){
      int layer = (*middleIndices)[i];
      //int octave = i / nOctaveLayers;
      vector< vector< int > > featureVectorParameters = fovea->getVectorsFeature();
      vector< int > bvec = featureVectorParameters[0];
      vector< int > etavec = featureVectorParameters[1];
      vector< int > levelvec = featureVectorParameters[2];
      int octave = etavec[layer];      
      if((*foveaLevel)[layer] == -1) continue;
      findMaximaInLayer( *sum, *mask_sum, *dets, *traces, *sizes,
			 *keypoints, octave, layer, hessianThreshold,
			 (*sampleSteps)[layer],
			 *fovea, (*margin)[layer], (*foveaLevel)[layer] );
    }
  }
  
  const Mat *sum;
  const Mat *mask_sum;
  const vector<Mat>* dets;
  const vector<Mat>* traces;
  const vector<int>* sizes;
  const vector<int>* sampleSteps;
  const vector<int>* middleIndices;
  vector<KeyPoint>* keypoints;
  int nOctaveLayers;
  float hessianThreshold;
  const vector<int> *foveaLevel;
  const vector<int> *margin;
  K* fovea;
  
  static Mutex findMaximaInLayer_m;
};

template <typename T, typename K>
Mutex SURFFindInvoker< T, K >::findMaximaInLayer_m;


/*
 * Find the maxima in the determinant of the Hessian in a layer of the
 * scale-space pyramid
 */
template <typename T, typename K>
void
SURFFindInvoker< T, K >::findMaximaInLayer( const Mat& sum, const Mat& mask_sum,
					    const vector<Mat>& dets, const vector<Mat>& traces,
					    const vector<int>& sizes, vector<KeyPoint>& keypoints,
					    int octave, int layer, float hessianThreshold, int sampleStep,
					    K fovea, int marginH, int foveaLevel ){
  cout << "find" << endl;
  // Wavelet Data
  const int NM=1;
  const int dm[NM][5] = { {0, 0, 9, 9, 1} };
  SurfHF Dm;
  
  int size = sizes[layer];
  
  // The integral image 'sum' is one pixel bigger than the source image
  int layer_rows = (sum.rows-1)/sampleStep;
  int layer_cols = (sum.cols-1)/sampleStep;
  
  // Ignore pixels without a 3x3x3 neighbourhood in the layer above
  int margin = (sizes[layer+1]/2)/sampleStep+1;
  
  if( !mask_sum.empty() )
    resizeHaarPattern( dm, &Dm, NM, 9, size, mask_sum.cols );
  
  int step = (int)(dets[layer].step/dets[layer].elemSize());
  
  for( int i = margin; i < layer_rows - margin; i++ ) {
    const float* det_ptr = dets[layer].ptr<float>(i);
    const float* trace_ptr = traces[layer].ptr<float>(i);
    for( int j = margin; j < layer_cols-margin; j++ ){
      float val0 = det_ptr[j];
      if( val0 > hessianThreshold ){
	//Coordinates for the start of the wavelet in the sum image. There
	//is some integer division involved, so don't try to simplify this
	//(cancel out sampleStep) without checking the result is the same
	int sum_i = sampleStep*(i-(size/2)/sampleStep);
	int sum_j = sampleStep*(j-(size/2)/sampleStep);
	
	//The 3x3x3 neighbouring samples around the maxima.
	//The maxima is included at N9[1][4]
	const float *det1 = &dets[layer-1].at<float>(i, j);
	const float *det2 = &dets[layer].at<float>(i, j);
	const float *det3 = &dets[layer+1].at<float>(i, j);
	float N9[3][9] = { { det1[-step-1], det1[-step], det1[-step+1],
			     det1[-1]  , det1[0] , det1[1],
			     det1[step-1] , det1[step] , det1[step+1]  },
			   { det2[-step-1], det2[-step], det2[-step+1],
			     det2[-1]  , det2[0] , det2[1],
			     det2[step-1] , det2[step] , det2[step+1]  },
			   { det3[-step-1], det3[-step], det3[-step+1],
			     det3[-1]  , det3[0] , det3[1],
			     det3[step-1] , det3[step] , det3[step+1]  } };
	
	//Check the mask - why not just check the mask at the center of the wavelet?
	if( !mask_sum.empty() ) {
	  const int* mask_ptr = &mask_sum.at<int>(sum_i, sum_j);
	  float mval = calcHaarPattern( mask_ptr, &Dm, 1 );
	  if( mval < 0.5 )
	    continue;
	}
	
	//Non-maxima suppression. val0 is at N9[1][4]
	if( val0 > N9[0][0] && val0 > N9[0][1] && val0 > N9[0][2] &&
	    val0 > N9[0][3] && val0 > N9[0][4] && val0 > N9[0][5] &&
	    val0 > N9[0][6] && val0 > N9[0][7] && val0 > N9[0][8] &&
	    val0 > N9[1][0] && val0 > N9[1][1] && val0 > N9[1][2] &&
	    val0 > N9[1][3]                    && val0 > N9[1][5] &&
	    val0 > N9[1][6] && val0 > N9[1][7] && val0 > N9[1][8] &&
	    val0 > N9[2][0] && val0 > N9[2][1] && val0 > N9[2][2] &&
	    val0 > N9[2][3] && val0 > N9[2][4] && val0 > N9[2][5] &&
	    val0 > N9[2][6] && val0 > N9[2][7] && val0 > N9[2][8] ) {
	  //Calculate the wavelet center coordinates for the maxima
	  float center_i = sum_i + (size-1)*0.5f;
	  float center_j = sum_j + (size-1)*0.5f;
	  
	  KeyPoint kpt( center_j, center_i, (float)sizes[layer],
			-1, val0, octave, (trace_ptr[j] > 0) - (trace_ptr[j] < 0) );
	  
	  // Interpolate maxima location within the 3x3x3 neighbourhood
	  int ds = size - sizes[layer-1];
	  int interp_ok = interpolateKeypoint( N9, sampleStep, sampleStep, ds, kpt );
	  
	  // Sometimes the interpolation step gives a negative size etc.
	  if( interp_ok  ){
	    //printf( "KeyPoint %f %f %d\n", point.pt.x, point.pt.y, point.size );
	    cv::AutoLock lock(findMaximaInLayer_m);
	    keypoints.push_back(kpt);
	  }
	}
      }
    }
  }
}

struct KeypointGreater
{
    inline bool operator()(const KeyPoint& kp1, const KeyPoint& kp2) const
    {
        if(kp1.response > kp2.response) return true;
        if(kp1.response < kp2.response) return false;
        if(kp1.size > kp2.size) return true;
        if(kp1.size < kp2.size) return false;
        if(kp1.octave > kp2.octave) return true;
        if(kp1.octave < kp2.octave) return false;
        if(kp1.pt.y < kp2.pt.y) return false;
        if(kp1.pt.y > kp2.pt.y) return true;
        return kp1.pt.x < kp2.pt.x;
    }
};

//////////////////
//// Updating ////
//////////////////


template <typename T, typename K>
static void fastFoveatedHessianDetector( const Mat& sum, const Mat& mask_sum, vector< KeyPoint >& keypoints, K fovea ){
  /* Sampling step along image x and y axes at first octave. This is doubled
     for each additional octave. WARNING: Increasing this improves speed,
     however keypoint extraction becomes unreliable. */
  const int SAMPLE_STEP0 = 1;

  vector< double > featureParameters = fovea.getParametersFeature();
  
  int nOctaveLayers = (int)featureParameters[0];
  float hessianThreshold = featureParameters[1];

  vector< vector< int > > featureVectorParameters = fovea.getVectorsFeature();
  vector< int > bvec = featureVectorParameters[0];
  vector< int > etavec = featureVectorParameters[1];
  vector< int > levelvec = featureVectorParameters[2];

  int nTotalLayers = (nOctaveLayers+2)*(bvec.size());
  int nMiddleLayers = nOctaveLayers*(bvec.size());

  vector<Mat> dets(nTotalLayers);
  vector<Mat> traces(nTotalLayers);
  vector<int> sizes(nTotalLayers);
  vector<int> sampleSteps(nTotalLayers);
  vector<int> middleIndices(nMiddleLayers);

  vector<int> foveaLevel(nTotalLayers);
  vector<int> margin(nTotalLayers);
  
  keypoints.clear();

  // Allocate space and calculate properties of each layer
  int index = 0, middleIndex = 0, step = SAMPLE_STEP0;

  vector< T > modelParameters = fovea.getParameters();

  for(unsigned int i = 0; i < (bvec.size()); i++) {
    for( int layer = 0; layer < nOctaveLayers+2; layer++ ) {
      // The integral image sum is one pixel bigger than the source image
      margin[index] = ((SURF_HAAR_SIZE0+SURF_HAAR_SIZE_INC*(nOctaveLayers+1))<<(etavec[i]-1))/2;
      if( bvec[i] == 0 )
	foveaLevel[index] = -1;
      else
	foveaLevel[index] = levelvec[i];
      dets[index].create( (modelParameters[2]).y, (modelParameters[2]).x, CV_32F );
      traces[index].create( (modelParameters[2]).y, (modelParameters[2]).x, CV_32F );
      sizes[index] = (SURF_HAAR_SIZE0 + SURF_HAAR_SIZE_INC*layer) << ( (featureVectorParameters[1])[i] - 1);
      sampleSteps[index] = 1 << (etavec[i] - 1);

      if( 0 < layer && layer <= nOctaveLayers )
	middleIndices[middleIndex++] = index;
      /*std::cout << index << " " << layer << ", sampleStep = " << sampleSteps[index] << "\t";
      std::cout << "Size: " << sizes[index] << ", eta = " << etavec[i] << std::endl;
      std::cout << "Size: " << sizes[index] << std::endl;
      std::cout << "Margin: " << margin[index] << std::endl; */
      index++;
    }
  }
  
  // Calculate hessian determinant and trace samples in each layer
  parallel_for_( Range(0, nTotalLayers),
		 SURFBuildInvoker< T, K >(sum, sizes, sampleSteps, dets, traces, fovea, margin, foveaLevel ) );
  
  // Find maxima in the determinant of the hessian
  parallel_for_( Range(0, nMiddleLayers),
		 SURFFindInvoker< T, K >(sum, mask_sum, dets, traces, sizes,
					 sampleSteps, middleIndices, keypoints,
					 nOctaveLayers, hessianThreshold, fovea, margin, foveaLevel ) );
  
  std::sort(keypoints.begin(), keypoints.end(), KeypointGreater());
}


//
// This function is used to apply Foveated Hessian Detector
//
template <typename T, typename K>
static void foveatedHessianDetector( InputArray _img, InputArray _mask, vector<KeyPoint>& keypoints, K fovea ) {
  Mat sum, mask1, msum;
  Mat img = _img.getMat();
  Mat mask = _mask.getMat();

  integral(img, sum, CV_32S);
  if(!mask.empty()) {
    cv::min(mask, 1, mask1);
    integral(mask1, msum, CV_32S);
  }

  CV_Assert(!img.empty() && img.depth() == CV_8U);
  if( img.channels() > 1 )
    cvtColor(img, img, COLOR_BGR2GRAY);

  vector< T > foveaParameters = fovea.getParameters();
  int m = foveaParameters[0].x;
  vector< Level< T > > levels;
  for ( int l = 0; l < m + 1; l++ )
    levels.push_back( fovea.getLevelFromFovea( l ) );
  
  fastFoveatedHessianDetector< T, K >(sum, msum, keypoints, fovea);
}



#endif
