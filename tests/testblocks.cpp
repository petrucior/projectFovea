#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../shape.hpp"
#include "../blocks.hpp"

using namespace cv;

int main(){
  // Testing constructor and print functions
  /*std::vector< Point > boundingBox;
  boundingBox.push_back( Point(0, 0) );
  boundingBox.push_back( Point(10, 10) );
  Blocks< Point > b( boundingBox );
  b.printVertices();
  std::vector< Point > boundingBoxVector;
  boundingBoxVector.push_back( Point(0, 0) );
  boundingBoxVector.push_back( Point(3, 3) );
  boundingBoxVector.push_back( Point(3, 3) );
  boundingBoxVector.push_back( Point(7, 7) );
  b.breakBlocks( boundingBoxVector );
  b.printVertices();*/
  
  std::vector< Point > boundingBoxA1, boundingBoxA2, boundingBoxA3, boundingBoxA4,
    boundingBoxA5, boundingBoxA6, boundingBoxA7, boundingBoxA8, boundingBoxB;
  // x eh linhas - contado de cima para baixo
  // y eh colunas - contado da esquerda para direita
  boundingBoxA1.push_back( Point(17, 17) );
  boundingBoxA1.push_back( Point(5, 5) );
  boundingBoxA2.push_back( Point(9, 17) );
  boundingBoxA2.push_back( Point(5, 5) );
  boundingBoxA3.push_back( Point(9, 9) );
  boundingBoxA3.push_back( Point(5, 5) );
  boundingBoxA4.push_back( Point(17, 9) );
  boundingBoxA4.push_back( Point(5, 5) );
  boundingBoxA5.push_back( Point(13, 17) );
  boundingBoxA5.push_back( Point(5, 5) );
  boundingBoxA6.push_back( Point(9, 13) );
  boundingBoxA6.push_back( Point(5, 5) );
  boundingBoxA7.push_back( Point(13, 9) );
  boundingBoxA7.push_back( Point(5, 5) );
  boundingBoxA8.push_back( Point(17, 13) );
  boundingBoxA8.push_back( Point(5, 5) );
  boundingBoxB.push_back( Point(13, 13) );
  boundingBoxB.push_back( Point(5, 5) );
  Blocks< Point > blockA1( boundingBoxA1 );
  Blocks< Point > blockA2( boundingBoxA2 );
  Blocks< Point > blockA3( boundingBoxA3 );
  Blocks< Point > blockA4( boundingBoxA4 );
  Blocks< Point > blockA5( boundingBoxA5 );
  Blocks< Point > blockA6( boundingBoxA6 );
  Blocks< Point > blockA7( boundingBoxA7 );
  Blocks< Point > blockA8( boundingBoxA8 );
  Blocks< Point > blockB( boundingBoxB );
  std::vector< Point > intersect;
  int direction;
  
  std::cout << "Analise entre A1 e B" << std::endl;
  blockA1.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA1, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[0].x << ", " << intersect[0].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A2 e B" << std::endl;
  blockA2.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA2, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[1].x << ", " << intersect[1].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A3 e B" << std::endl;
  blockA3.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA3, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[2].x << ", " << intersect[2].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  
  std::cout << "Analise entre A4 e B" << std::endl;
  blockA4.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA4, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[3].x << ", " << intersect[3].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A5 e B" << std::endl;
  blockA5.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA5, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[4].x << ", " << intersect[4].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A6 e B" << std::endl;
  blockA6.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA6, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[5].x << ", " << intersect[5].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A7 e B" << std::endl;
  blockA7.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA7, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[6].x << ", " << intersect[6].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::cout << "Analise entre A8 e B" << std::endl;
  blockA8.printVertices();
  blockB.printVertices();
  std::cout << blockB.intersectionShapeByVertices( blockA8, intersect, direction ) << std::endl;
  if ( intersect.size() == 0 ){
    std::cout << "without intersection " << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }
  else{
    std::cout << "point: " << intersect[7].x << ", " << intersect[7].y << std::endl;
    std::cout << "direction: " << direction << std::endl;
  }

  std::vector< Shape< Point >* > shapes;
  //shapes.push_back( &blockA1 );
  shapes.push_back( &blockA2 );
  //shapes.push_back( &blockA3 );
  shapes.push_back( &blockA4 );
  /*shapes.push_back( &blockA5 );
  shapes.push_back( &blockA6 );
  shapes.push_back( &blockA7 );
  shapes.push_back( &blockA8 );*/
  std::cout << " -------------" << std::endl;
  std::cout << "  breakBlocks " << std::endl; 
  std::cout << " -------------" << std::endl;
  blockB.breakBlocks( shapes );
  
  // Testing distance function ( It is necessary to add vertices in clockwise )
  /*Point vA = Point( 0, 0 );
  Point vB = Point( 10, 0 );
  Point p = Point( 5, 5 );
  std::cout << r.distance( vA, vB, p ) << std::endl;*/
  
  // Testing intersectionShape with point
  /*Point point = Point(3, 2);
  std::cout << r.intersectionShapeByPoint( point ) << std::endl;

  // Testing intersectionShape with another shape
  std::vector< Point > boundingBox2;
  std::vector< Point > points;
  boundingBox2.push_back( Point(0, 0) );
  boundingBox2.push_back( Point(4, 4) );
  Rectangle< Point > s( boundingBox2 );
  std::cout << r.intersectionShape( s, points ) << std::endl;*/
  return 0;
}
