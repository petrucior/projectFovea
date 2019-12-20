/*#include <iostream>
#include <string>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../statistics.hpp"*/

#include <iostream>
#include <string>
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include "../statistics.hpp"

using namespace cv;

static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  if ( event == EVENT_MOUSEMOVE ){
    Fovea< Point > *fovea = (Fovea< Point > *) _fovea;
    fovea->updateFovea( Point(x, y ) );
    //std::cout << "( " << x << ", " << y << " )" << std::endl;
  }
}

int main( int argc, char** argv ){
  
  // Testing constructor and print functions
  std::cout << argv[1] << std::endl;
  cv::Mat model = cv::imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  cv::Mat scene = cv::imread(argv[2]);//, CV_LOAD_IMAGE_GRAYSCALE):
  cv::Mat sceneFoveated;
  cv::Point fs = cv::Point( -60, 0 );
  cv::Scalar color = cv::Scalar( 255, 255, 255 );
  // Extraction features of model
  cv::Ptr<cv::FeatureDetector> detector;
  cv::Ptr<cv::DescriptorExtractor> descriptor;

  // -----------------------
  // Configuration to _SURF_
  // -----------------------
  detector = cv::xfeatures2d::SURF::create();
  descriptor = cv::xfeatures2d::SURF::create();
  Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(120, 120), cv::Point( scene.cols, scene.rows ), fs );
  
  // -----------------------
  // Configuration to _ORB_
  // -----------------------
  //detector = cv::ORB::create();
  //descriptor = cv::ORB::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(120, 120), cv::Point( scene.cols, scene.rows ), fs );
  
  // -----------------------
  // Configuration to _KAZE_
  // -----------------------
  //detector = cv::KAZE::create();
  //descriptor = cv::KAZE::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 4, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
    
  // -----------------------
  // Configuration to _AKAZE_
  // -----------------------
  //detector = cv::AKAZE::create();
  //descriptor = cv::AKAZE::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
  
  // -----------------------
  // Configuration to _BRISK_
  // -----------------------
  //detector = cv::BRISK::create();
  //descriptor = cv::BRISK::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
    
  
  // -----------------------
  // Configuration to Model
  // -----------------------
  std::vector< cv::KeyPoint > modelKeypoints;
  cv::Mat modelDescriptors;
  detector->detect ( model, modelKeypoints );
  descriptor->compute ( model, modelKeypoints, modelDescriptors );
    
  
  //Statistics< cv::Point2f > *s = new Statistics< cv::Point2f >();
  Statistics< cv::Point > *s = new Statistics< cv::Point >();
  //s->plotProportion( fscene, scene, model, _SURF_, modelKeypoints, modelDescriptors, 0.0, 1000.0 );
  
  cv::namedWindow("sceneFoveated", 1);
  while ( true ){
    setMouseCallback( "sceneFoveated", on_mouse, fscene );
    sceneFoveated = fscene->foveatedImage( scene, cv::Scalar(255, 255, 255) );
    imshow("sceneFoveated", sceneFoveated );
    drawKeypoints( model, modelKeypoints, model, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    imshow("modelFoveated", model );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'e' ){
      //fscene->foveatedFeatures( scene, _KAZE_, MRMF );
      //fscene->foveatedFeatures( scene, _ORB_, MRMF );
      //fscene->foveatedFeatures( scene, _BRISK_, MRMF );
      fscene->foveatedFeatures( scene, _SURF_, MRMF );
      fscene->matching( scene, model, modelKeypoints, modelDescriptors );
      
      std::cout << s->functionFovea( fscene, 15.0, 1000.0 ) << std::endl;
      
    }
  }
  
  /*double referencePotential = 0.8;
  std::vector< double > potentials;  
  potentials.push_back( 0.5 );
  potentials.push_back( 0.7 );
  potentials.push_back( 0.9 );
  potentials.push_back( 0.3 );
  int index = s->localGradient( referencePotential, potentials );
  cv::Point px = cv::Point( 10, 10 );
  std::vector< cv::Point2f > regionSearched, updatedRegion;
  regionSearched.push_back( cv::Point2f( 0, 0 ) );
  regionSearched.push_back( cv::Point2f( 512, 512 ) );
  int configuration = 0;
  updatedRegion = s->reduceRegionByLocalGradient( index, px, regionSearched, configuration );
  std::cout << std::endl;
  for( int i = 0; i < regionSearched.size(); i++ )
    std::cout << "( " << updatedRegion[i].x << ", " << updatedRegion[i].y << " ) " << std::endl;
  */

  /*std::vector< cv::Point2f > pointA, pointB;
  pointA.push_back( cv::Point2f( 1, 0 ) );
  //pointA.push_back( cv::Point2f( -1, -1 ) );
  pointA.push_back( cv::Point2f( 0, 0 ) );
  //pointB.push_back( cv::Point2f( 10, 0 ) );
  //pointB.push_back( cv::Point2f( 8, 2 ) );
  pointB.push_back( cv::Point2f( 9, 1 ) );
  pointB.push_back( cv::Point2f( 11, -1 ) );
  cv::Point2f pointEstimated = s->intersectionLocalGradient(pointA, pointB);
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;
  */
  
  /*std::vector< cv::Point2f > points;
  points.push_back( cv::Point2f( 4, 8 ) );
  points.push_back( cv::Point2f( 9, 6 ) );
  points.push_back( cv::Point2f( 5, 5 ) );
  points.push_back( cv::Point2f( 0, 0 ) );
  points.push_back( cv::Point2f( 7, 7 ) );
  // point to be estimated is (6, 6)
  std::vector< double > potentials;
  
  //      pos ----- potential  \  complement ----- radius
  // ref 6, 6 -----    1       \   0.0001    ----- 0.001
  //     3, 3 -----    0.5     \   0.5       ----- 5
  //     4, 8 -----    0.718   \   0.282     ----- 2.82
  //     9, 6 -----    0.7     \   0.3       ----- 3
  //     5, 5 -----    0.859   \   0.141     ----- 1.41
  //     0, 0 -----    0       \   1         ----- 8.48
  //     7, 7 -----    0.859   \   0.141     ----- 1.41
  
  potentials.push_back(0.718);
  potentials.push_back(0.7);
  potentials.push_back(0.859);
  potentials.push_back(0);
  potentials.push_back(0.859);
  cv::Point2f pointEstimated = s->maximumLikelihoodEstimator( points, potentials, 1 );
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;*/

  /*std::vector< cv::Point2f > points;
  points.push_back( cv::Point2f( 4, 8 ) );
  points.push_back( cv::Point2f( 9, 6 ) );
  points.push_back( cv::Point2f( 5, 5 ) );
  //points.push_back( cv::Point2f( 0, 0 ) );
  //points.push_back( cv::Point2f( 7, 7 ) );
  // point to be estimated is (6, 6)
  std::vector< double > invpotentials;

  //      pos ----- potential  \  complement ----- radius
  // ref 6, 6 -----    1       \   0.0001    ----- 0.001
  //     3, 3 -----    0.5     \   0.5       ----- 5
  //     4, 8 -----    0.718   \   0.282     ----- 2.82
  //     9, 6 -----    0.7     \   0.3       ----- 3
  //     5, 5 -----    0.859   \   0.141     ----- 1.41
  //     0, 0 -----    0       \   1         ----- 8.48
  //     7, 7 -----    0.859   \   0.141     ----- 1.41
  
  invpotentials.push_back(2.82);
  invpotentials.push_back(3);
  invpotentials.push_back(1.41);
  //invpotentials.push_back(8.48);
  //invpotentials.push_back(1.41);
  cv::Point2f pointEstimated = s->trilaterationEstimator(points, invpotentials);
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;*/
  
  
  return 0;
  
}
