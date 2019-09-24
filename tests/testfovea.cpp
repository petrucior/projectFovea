#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
#include "../fovea.hpp"

using namespace cv;

int main( int argc, char** argv ){
  // Testing constructor and print functions
  std::cout << argv[1] << std::endl;
  Mat img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  //Fovea< Point > f( img, 4, Point(30, 30), Point( 0, 0 ) );
  Fovea< Point > f( 4, Point(40, 40), Point( img.cols, img.rows ), Point( 0, 0 ) );
  f.foveatedFeatures( img, _KAZE_, MRMF );
  //img = l.get( img );
  //imshow( "image", img );
  //waitKey( 0 );
  return 0;
}
