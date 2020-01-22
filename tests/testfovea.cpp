#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
#include "../fovea.hpp"

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
  Mat img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  //Fovea< Point > f( img, 4, Point(30, 30), Point( 0, 0 ) );
  // Configuration to _KAZE_
  //Fovea< Point > *f = new Fovea< Point >( 4, Point(60, 60), Point( img.cols, img.rows ), Point( 0, 0 ) );
  // Configuration to _ORB_
  Fovea< Point > *f = new Fovea< Point >( 4, Point(80, 80), Point( img.cols, img.rows ), Point( 0, 0 ), _BLOCKS_ );
  //f.foveatedFeatures( img, _ORB_, MRMF );
  //f->foveatedFeatures( img, _KAZE_, MRMF );
  //img = l.get( img );
  namedWindow("image", 1);
  Mat foveated;
  while ( true ){
    setMouseCallback( "image", on_mouse, f );
    foveated = f->foveatedImage( img, cv::Scalar(255, 255, 255) );
    f->foveatedFeatures( img, _ORB_, MRMF, *f );
    //f->foveatedFeatures( img, _KAZE_, MRMF );
    imshow("image", foveated );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'u' ){
      f->foveatedFeatures( img, _KAZE_, MRMF, *f );
    }
  }
    

  return 0;
}
