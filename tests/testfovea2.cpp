#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
#include "../fovea.hpp"

using namespace cv;

static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  Fovea< Point > *fovea = (Fovea< Point > *) _fovea;
  fovea->updateFovea( Point(x, y ));
  std::cout << "atualiza" << std::endl;
}

int main( int argc, char** argv ){
  
  // Testing constructor and print functions
  std::cout << argv[1] << std::endl;
  Mat img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  //Fovea< Point > f( img, 4, Point(30, 30), Point( 0, 0 ) );
  Fovea< Point > *f = new Fovea< Point >( 4, Point(100, 100), Point( img.cols, img.rows ), Point( 0, 0 ) );
  //f.foveatedFeatures( img, _ORB_, MRMF );
  f->foveatedFeatures( img, _KAZE_, MRMF );
  //img = l.get( img );
  //imshow( "image", img );
  //waitKey( 0 );
    

  return 0;
}
