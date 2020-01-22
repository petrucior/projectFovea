#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
//#include "../level.hpp"
//#include "../fovea.hpp"
#include "../multifovea.hpp"

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
  Mat img = imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  Mat multifoveated;
  std::vector< Point > fs;
  // Foveas centralized position
  /*fs.push_back( Point( 0, 0 ) );
  fs.push_back( Point( 0, 0 ) );
  fs.push_back( Point( 0, 0 ) );
  fs.push_back( Point( 0, 0 ) );
  fs.push_back( Point( 0, 0 ) );*/
  
  // Foveas simetrically positioned
  fs.push_back( Point( -40, 40 ) );
  fs.push_back( Point( 40, 40 ) );
  fs.push_back( Point( 40, -40 ) );
  fs.push_back( Point( -40, -40 ) );
  std::vector<cv::Scalar> colors;
  srand (time(NULL)); // Initialize random seed                                                                                                                                                 
  for (int i = 0; i < fs.size(); i++){
    int r = rand() % 256; // 0 - 255                                                                                                                                                            
    int g = rand() % 256; // 0 - 255                                                                                                                                                            
    int b = rand() % 256; // 0 - 255                                                                                                                                                            
    colors.push_back(cv::Scalar(b, g, r));
  }
  // Configuration to _KAZE_
  Multifovea< Point > *m = new Multifovea< Point >( 4, Point(20, 20), Point( img.cols, img.rows ), fs, _REEXECUTION_, _BLOCKS_ );
  while ( true ){
    //setMouseCallback( "image", on_mouse, m->getFovea(0) );
    multifoveated = m->multifoveatedImage( img );
    //multifoveated = m->multifoveaLevelsImage( img, colors );
    imshow("image", multifoveated );
    waitKey( 0 );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'u' )
      m->foveatedFeatures( img, _KAZE_, MRMF );
  }
  
  return 0;
}
