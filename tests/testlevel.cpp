#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "../shape.hpp"
//#include "../rectangle.hpp"
#include "../level.hpp"

using namespace cv;

int main( int argc, char** argv ){
  // Testing constructor and print functions
  /*std::vector< Point > boundingBox;
  boundingBox.push_back( Point(0, 0) );
  boundingBox.push_back( Point(10, 10) );
  //Rectangle< Point > *r = new Rectangle< Point>( boundingBox );
  //Shape< Point >* s = dynamic_cast< Shape< Point >* >( r );
  //r->printVertices();
  //
  // Two options between constructors
  //
  std::cout << argv[1] << std::endl;
  Mat img = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  //Level< Point > l( 0, r );
  Level< Point > l( 0, 4, Point(30, 30), Point( img.cols, img.rows ), Point( 0, 0 ) ); 
  //boundingBox = l.boundingBox( 4, 4, Point(30, 30), Point( 100, 100 ), Point( 0, 0 ) );
  //std::cout << boundingBox[0] << ", " << boundingBox[1] << std::endl;
  img = l.getLevel( img );
  imshow( "image", img );
  waitKey( 0 );*/

  std::cout << argv[1] << std::endl;
  Mat img = imread( argv[1], CV_LOAD_IMAGE_GRAYSCALE );
  Level< Point > A1( 4, 4, Point(100, 100), Point( img.cols, img.rows ), Point( -40, 40 ), _BLOCKS_ );
  Level< Point > A2( 4, 4, Point(100, 100), Point( img.cols, img.rows ), Point( 40, 40 ), _BLOCKS_ );
  Level< Point > A3( 4, 4, Point(100, 100), Point( img.cols, img.rows ), Point( 40, -40 ), _BLOCKS_ );
  Level< Point > A4( 4, 4, Point(100, 100), Point( img.cols, img.rows ), Point( -40, -40 ), _BLOCKS_ );
  Level< Point > B( 0, 4, Point(100, 100), Point( img.cols, img.rows ), Point( 0, 0 ), _BLOCKS_ );

  /*Level< Point > l( 0, 4, Point(100, 100), Point( img.cols, img.rows ), Point( 0, 0 ), _BLOCKS_ );
  img = l.getLevel( img );
  imshow( "image", img );
  waitKey( 0 );*/
  /*vector< Mat > vecMat = l.getLevel2( img );
  for ( int i = 0; i < vecMat.size(); i++ ){
    imshow( "image", vecMat[i] );
    waitKey( 0 );
  }*/

  imshow("image real", img );

  Mat imgref = B.getLevel( img );
  imshow( "image reference", imgref );
    
  vector< Level< Point > > levels;
  levels.push_back( A1 );
  //levels.push_back( A2 );
  //levels.push_back( A3 );
  //levels.push_back( A4 );
  B.updateLevel( levels );

  vector< Mat > vecMat = B.getLevel2Render( img );
  cout << "vecMat: " << vecMat.size() << endl;
  for ( int i = 0; i < vecMat.size(); i++ ){
    std::cout << "image " << i << std::endl;
    imshow( "image", vecMat[i] );
    waitKey( 0 );
  }
  
  return 0;
}
