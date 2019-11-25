#include <iostream>
#include <string>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
#include "../fovea.hpp"
//#include "../multifovea.hpp"

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
  Mat model = imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  Mat scene = imread(argv[2]);//, CV_LOAD_IMAGE_GRAYSCALE):
  namedWindow("sceneFoveated", 1);
  Mat sceneFoveated;
  Point fs = Point( -60, 0 );
  cv::Scalar color = cv::Scalar( 255, 255, 255 );
  // Extraction features of model
  cv::Ptr<cv::FeatureDetector> detector;
  cv::Ptr<cv::DescriptorExtractor> descriptor;
  // Configuration to _ORB_
  detector = cv::ORB::create();
  descriptor = cv::ORB::create();
  
  Fovea< Point > *fscene = new Fovea< Point >( 2, Point(120, 120), Point( scene.cols, scene.rows ), fs );

  // Configuration to _KAZE_
  //detector = cv::KAZE::create();
  //descriptor = cv::KAZE::create();
  std::vector< cv::KeyPoint > modelKeypoints;
  cv::Mat modelDescriptors;
  detector->detect ( model, modelKeypoints );
  descriptor->compute ( model, modelKeypoints, modelDescriptors );
  
  //Fovea< Point > *fscene = new Fovea< Point >( 4, Point(100, 100), Point( scene.cols, scene.rows ), fs );
  
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
      fscene->foveatedFeatures( scene, _ORB_, MRMF );
      fscene->matching( scene, model, modelKeypoints, modelDescriptors );
      /*char buffer [33];
      double ratio = fscene->matching2( scene, model, modelKeypoints, modelDescriptors );
      std::cout << ratio << std::endl;*/
    }
  }
  
  return 0;
}
