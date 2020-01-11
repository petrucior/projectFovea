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

using namespace cv;

int main( int argc, char** argv ){
  
  // Testing constructor and print functions
  Mat img = imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  String ymlFile = argv[2];
  imshow( "image", img );
  waitKey( 0 );
  
  return 0;
}
