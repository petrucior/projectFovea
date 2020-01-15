/*#include <iostream>AOA
#include <string>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../statistics.hpp"*/

#include <iostream>
#include <stdlib.h>
#include <string>
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

#include "../fovea.hpp"
//#include "../foveatedHessianDetector.hpp"

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  if ( event == EVENT_MOUSEMOVE ){
    Fovea< Point2f > *fovea = (Fovea< Point2f > *) _fovea;
    fovea->updateFovea( Point2f(x, y ) );
    //cout << "( " << x << ", " << y << " )" << endl;
  }
}

int main( int argc, char** argv ){
  Mat model = imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  Mat scene = imread(argv[2]);//, CV_LOAD_IMAGE_GRAYSCALE):
  String ymlFile = argv[3];
  Mat sceneFoveated = scene;
  Mat multifoveated = scene;
  vector< Point2f > fs;
  // Foveas simetrically positioned
  fs.push_back( Point2f( -40.0, 40.0 ) );
  vector< Scalar > colors;
  srand (time(NULL)); // Initialize random seed                                                                                                                                                 
  for (int i = 0; i < fs.size(); i++){
    int r = rand() % 256; // 0 - 255                                                                                                                                                            
    int g = rand() % 256; // 0 - 255                                                                                                                                                            
    int b = rand() % 256; // 0 - 255                                                                                                                                                            
    colors.push_back( Scalar(b, g, r) );
  }
  // Extraction features of model
  Ptr< FeatureDetector > detector;
  Ptr< DescriptorExtractor > descriptor;

  // -----------------------
  // Configuration to _SURF_
  // -----------------------
  detector = SURF::create();
  descriptor = SURF::create();
  //Fovea< Point2f > *fscene = new Fovea< Point2f >( 4, Point2f(120, 120), Point2f( scene.cols, scene.rows ), fs[0] );
  //Fovea< Point2f > fscene( 4, Point2f(120, 120), Point2f( scene.cols, scene.rows ), fs[0] );
  Fovea< Point2f > fscene( scene, ymlFile, 0);
  
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
  
  
  cv::namedWindow("sceneFoveated", 1);
  while ( true ){
    setMouseCallback( "sceneFoveated", on_mouse, &fscene );
    //sceneFoveated = fscene.foveatedImage( scene, colors[0], MMF );
    imshow("sceneFoveated", sceneFoveated );
    drawKeypoints( model, modelKeypoints, model, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    imshow("modelFoveated", model );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'e' ){
      //fscene.foveatedFeatures( scene, _SURF_, MRMF, fscene );
      fscene.foveatedFeatures( scene, _SURF_, MMF, fscene );
      //fscene->matching( scene, model, modelKeypoints, modelDescriptors );
      sceneFoveated = fscene.foveatedImage( scene, colors[0], MMF );
    }
  }
  
  return 0;
  
}
