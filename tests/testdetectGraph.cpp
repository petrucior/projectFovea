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
  
  FILE *file;
  char buffer[100];
  file = fopen ("data/graph.dat","w");
  if ( file != NULL ){
    fprintf( file, "#%c #%c #%c \n", 'x', 'y', 'z');
    for ( int i = 0; i < scene.rows; i+=10 ){
      for ( int j = 0; j < scene.cols; j+=10 ){
	fscene->updateFovea( Point( i, j ) );
	fscene->foveatedFeatures( scene, _ORB_, MRMF );
	double ratio = fscene->matching2( scene, model, modelKeypoints, modelDescriptors );
	//std::cout << ratio << std::endl;
	fprintf( file, "%d %d %f \n", i, j, ratio );
      }
    }
  }
  fclose( file );
  
  return 0;
}
