#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
//#include "../fovea.hpp"
#include "../multifovea.hpp"

using namespace cv;

/*static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  if ( event == EVENT_MOUSEMOVE ){
    Fovea< Point > *fovea = (Fovea< Point > *) _fovea;
    fovea->updateFovea( Point(x, y ) );
    //std::cout << "( " << x << ", " << y << " )" << std::endl;
  }
}*/

int main( int argc, char** argv ){
  
  // Testing constructor and print functions
  std::cout << argv[1] << std::endl;
  Mat img = imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  namedWindow("image", 1);
  Mat multifoveated;
  std::vector< Point > fs;
  fs.push_back( Point( -40, 40 ) );
  fs.push_back( Point( 40, 40 ) );
  fs.push_back( Point( 40, -40 ) );
  fs.push_back( Point( -40, -40 ) );
  Multifovea< Point > *m = new Multifovea< Point >( 10, Point(20, 20), Point( img.cols, img.rows ), fs, REEXECUTION );
  while ( true ){
    //setMouseCallback( "image", on_mouse, f );
    multifoveated = m->multifoveatedImage( img );
    imshow("image", multifoveated );
    waitKey( 0 );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
  }
  
  return 0;
}
