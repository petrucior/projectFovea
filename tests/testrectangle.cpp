#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../shape.hpp"
#include "../rectangle.hpp"

using namespace cv;

int main(){
  // Testing constructor and print functions
  std::vector< Point > boundingBox;
  boundingBox.push_back( Point(0, 0) );
  boundingBox.push_back( Point(10, 10) );
  Rectangle< Point > r( boundingBox );
  r.printVertices();
  
  // Testing distance function ( It is necessary to add vertices in clockwise )
  /*Point vA = Point( 0, 0 );
  Point vB = Point( 10, 0 );
  Point p = Point( 5, 5 );
  std::cout << r.distance( vA, vB, p ) << std::endl;*/
  
  // Testing intersectionShape with point
  Point point = Point(3, 2);
  std::cout << r.intersectionShapeByPoint( point ) << std::endl;

  // Testing intersectionShape with another shape
  std::vector< Point > boundingBox2;
  std::vector< Point > points;
  boundingBox2.push_back( Point(0, 0) );
  boundingBox2.push_back( Point(4, 4) );
  Rectangle< Point > s( boundingBox2 );
  int direction;
  std::cout << r.intersectionShape( s, points, direction ) << std::endl;
  return 0;
}
