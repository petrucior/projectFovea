/**
 * \file blocks.hpp
 *
 * \brief This file contains the prototype of shape blocks.
 *
 * \author 
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date Jan 2020
 *
 * This file is part of projectFovea software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License, or (at your option) any later
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BLOCKS_HPP
#define BLOCKS_HPP

#include <iostream> //std::cout, std::endl
#include <algorithm> // std::max, std::min
//#include <stdlib.h> // abs
//#include "shape.hpp" //Shape
#include <vector> //std::vector
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

using namespace std;
using namespace cv;

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

enum directions {SOUTHEAST, NORTHEAST, SOUTHWEST, NORTHWEST};

/**
 * \class Blocks
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea.
 *
 * \tparam T - Generic representation for type Point
 */
template < typename T > // cv::Point 
class Blocks : public Shape< T > {
public:
  //
  // Methods
  //
  /**
   * \fn Blocks( vector< T > boundingBox )
   *
   * \brief Constructor default.
   * This method initialize the start and final vertices of the shape.
   *
   * \param boundingBox - Vector containing 2 positions with tuple
   * information the limits of rectangular region ( delta and size )
   * to create the fovea
   */
  Blocks( vector< T > boundingBox );

  /**
   * \fn void breakBlocks( vector< Shape< T >* > shapes );
   *
   * \brief Constructor breaks the lowerBound/upperBound of current shape 
   * in blocks and update the vertices of shape.
   *
   * \param shapes - Vector contains n shapes ( Note: All shapes are on 
   * the same level )
   */
  void breakBlocks( vector< Shape< T >* > shapes );
  
  /**
   * \fn ~Blocks()
   *
   * \brief Destructor default
   */
  ~Blocks();
  
  /**
   * \fn virtual void printVertices()
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method print all vertices of the shape.
   */
  void printVertices();
    
  /**
   * \fn bool intersectionShapeByVertices( const Shape< T >& shape, const vector< T >& pointIntersection, int& direction )
   * 
   * \brief This method checks if the shape intersect the region delimited  by vertices and return
   * the point of intersection. Consider x ( lines - from top to bottom ) and y ( columns - from 
   * left to right )
   *
   * \param shape - Shape analyzed
   * \param pointIntersection - reference vertex of the another block, assumption: they intersect
   * \param direction - 0 (southeast), 1 (northeast), 2 (southwest), 3 (northwest)
   *
   * \return True if shape intersect and false otherwise.
   */
  bool intersectionShapeByVertices( Shape< T >& shape, vector< T >& pointIntersection, int& direction );

  /**
   * \fn bool intersectionShape( const Shape< T >& shape, RefVertex< T >& v )
   * 
   * \brief This method checks if the shape intersect the region delimited  by vertices and return
   * the point of intersection. Consider x ( lines - from top to bottom ) and y ( columns - from 
   * left to right )
   *
   * \param shape - Shape analyzed
   * \param v - Reference to vertex of intersection and its direction
   *
   * \return True if shape intersect and false otherwise.
   */
  bool intersectionShape( Shape< T >& shape, RefVertex< T >& v );
  
private:
  //
  // Attributes
  //
  T delta;
  T size;

  /**
   * \fn bool compareVertex( RefVertex< T > va, RefVertex< T > vb )
   * 
   * \brief This method compares two references vertex to be used in
   * the sort algorithm
   *
   * \param va, vb - References vertex that contains intersectionPoint 
   * and direction associated
   *
   *\return True if the x-axis of va is smaller than x-axis of vb and 
   * false otherwise.
   */
  static bool compareVertex( RefVertex< T > va, RefVertex< T > vb );

};

#endif

/**
 * \fn Blocks( vector< T > boundingBox )
 *
 * \brief Constructor default.
 * This method initialize the start and final vertices of the shape.
 *
 * \param boundingBox - Vector containing 2 positions with tuple
 * information the limits of rectangular region ( delta and size )
 * to create the fovea
 */
template <typename T>
Blocks< T >::Blocks( vector< T > boundingBox ) : Shape< T >( boundingBox ){
  delta = boundingBox[0];
  size = boundingBox[1];
  this->updateVertices( boundingBox );
  this->boundingBoxShape = boundingBox;
}

/**
 * \fn void breakBlocks( vector< Shape< T >* > shapes );
 *
 * \brief Constructor breaks the lowerBound/upperBound of current shape 
 * in blocks and update the vertices of shape.
 *
 * \param shapes - Vector contains n shapes ( Note: All shapes are on 
 * the same level )
 */
template <typename T>
void
Blocks< T >::breakBlocks( vector< Shape< T >* > shapes ){
  vector< RefVertex< T > > vertexVector;
  for( int s = 0; s < shapes.size(); s++ ){
    RefVertex< T > vertex;
    if ( intersectionShape( (dynamic_cast<Blocks &>(*shapes[s])), vertex ) )
      vertexVector.push_back( vertex );
  }
  if ( vertexVector.size() == 0 ) return;
  else{
    //creates (t+1) regions
    sort( vertexVector.begin(), vertexVector.end(), compareVertex );
    vector<int> lowerBound;
    vector<int> upperBound;
    for( int bound = 0; bound <= vertexVector.size(); bound++ ) {
      lowerBound.push_back(delta.y);
      upperBound.push_back(delta.y+size.y);
    }
#ifdef DEBUG
    for( int bound = 0; bound <= vertexVector.size(); bound++ ) {
      printf("block %d = [%d %d]\n", bound, lowerBound[bound], upperBound[bound]);
    }
#endif
    //adjusts the lower and upperBound for each region
    for( int i = 0; i < vertexVector.size(); i++ ) {
      if(vertexVector[i].direction == NORTHWEST) {
	for(int j = 0; j <= i; j++)
	  lowerBound[j] = max(lowerBound[j], vertexVector[i].intersectionPoint.y);
      } else if(vertexVector[i].direction == NORTHEAST) {
	for(int j = i+1; j <= vertexVector.size(); j++)
	  lowerBound[j] = max(lowerBound[j], vertexVector[i].intersectionPoint.x);
      } else if(vertexVector[i].direction == SOUTHEAST) {
	for(int j = i+1; j <= vertexVector.size(); j++)
	  upperBound[j] = min(upperBound[j], vertexVector[i].intersectionPoint.y);
      } else {
	for(int j = 0; j <= i; j++) // SOUTHWEST
	  upperBound[j] = min(upperBound[j], vertexVector[i].intersectionPoint.x);
      }
    }
#ifdef DEBUG
    for( int bound = 0; bound <= vertexVector.size(); bound++ ) {
      printf("block updated %d = [%d %d]\n", bound, lowerBound[bound], upperBound[bound]);
    }
#endif
    
    vector< T > newBoundingBox;
    T newDelta = T( lowerBound[0], delta.y );
    T newSize = T( upperBound[0] - lowerBound[0], vertexVector[0].intersectionPoint.y - newDelta.y );  
    if ( ( newSize.x > 0 ) && ( newSize.y > 0 ) ){
      newBoundingBox.push_back( newDelta );
      newBoundingBox.push_back( newSize );
    }
    for(int i = 0; i < vertexVector.size(); i++) {
      newDelta = T( lowerBound[i+1], vertexVector[i].intersectionPoint.y );
      if ( i == vertexVector.size() - 1 )
	newSize.y = (delta.y + size.y) - newDelta.y;
      else
	newSize.y = vertexVector[i+1].intersectionPoint.y - newDelta.y;
      newSize.x = upperBound[i+1] - lowerBound[i+1];
      if ( ( newSize.x > 0 ) && ( newSize.y > 0 ) ){
	newBoundingBox.push_back( newDelta );
	newBoundingBox.push_back( newSize );
      }
    }
    this->updateVertices( newBoundingBox );
  }  
}

/**
 * \fn ~Blocks()
 *
 * \brief Destructor default
 */
template <typename T>
Blocks< T >::~Blocks(){
  vector< T >().swap( this->vertices ); // Free the memoyr
  vector< T >().swap( this->boundingBoxShape ); // Free the memory
}

/**
 * \fn void printVertices()
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * This method print all vertices of the shape.
 */
template <typename T>
void 
Blocks< T >::printVertices(){
  //cout << this->vertices.size() << endl;
  for ( int v = 0; v < this->vertices.size(); v+=2 )
    cout << "delta: ( " << this->vertices[v].x <<  ", " << this->vertices[v].y << "), size: (" << this->vertices[v+1].x << ", " << this->vertices[v+1].y  << ")"<< endl;
}

/**
 * \fn bool intersectionShapeByVertices( const Shape< T >& shape, const vector< T >& pointIntersection, int& direction )
 *
 * \brief This method checks if the shape intersect the region delimited  by vertices and return
 * the point of intersection. Consider x ( lines - from top to bottom ) and y ( columns - from 
 * left to right )
 *
 * \param shape - Shape analyzed
 * \param pointIntersection - reference vertex of the another block, assumption: they intersect
 * \param direction - 0 (southeast), 1 (northeast), 2 (southwest), 3 (northwest)
 *
 * \return True if shape intersect and false otherwise.
 */
template <typename T>
bool 
Blocks< T >::intersectionShapeByVertices( Shape< T >& shape, vector< T >& pointIntersection, int& direction ){
  vector< T > boundingBoxShape = shape.getBoundingBox();
  T deltaShape = boundingBoxShape[0];
  T sizeShape = boundingBoxShape[1];
  int intersectx = - max( delta.x, deltaShape.x ) + min( delta.x + size.x, deltaShape.x + sizeShape.x );
  int intersecty = - max( delta.y, deltaShape.y ) + min( delta.y + size.y, deltaShape.y + sizeShape.y );
  T point;
  if ( ( intersectx > 0 ) || ( intersecty > 0 ) ){
    direction = 0; // its southeast at first
    if ( delta.x > deltaShape.x ){
      point.x = deltaShape.x + sizeShape.x;
      direction += 1; // direction is changed to north
    }
    else
      point.x = deltaShape.x;
    if ( delta.y > deltaShape.y ){
      point.y = deltaShape.y + sizeShape.y;
      direction += 2; // direction is changed to west
    }
    else
      point.y = deltaShape.y;
    
    pointIntersection.push_back( point );
  }
  if ( ( intersectx > 0 ) || ( intersecty > 0 ) ) return true;
  return false;
}

/**
 * \fn bool intersectionShape( const Shape< T >& shape, RefVertex< T >& v )
 * 
 * \brief This method checks if the shape intersect the region delimited  by vertices and return
 * the point of intersection. Consider x ( lines - from top to bottom ) and y ( columns - from 
 * left to right )
 *
 * \param shape - Shape analyzed
 * \param v - Reference to vertex of intersection and its direction
 *
 * \return True if shape intersect and false otherwise.
 */
template <typename T>
bool
Blocks< T >::intersectionShape( Shape< T >& shape, RefVertex< T >& v ){
  vector< T > boundingBoxShape = shape.getBoundingBox();
  T deltaShape = boundingBoxShape[0];
  T sizeShape = boundingBoxShape[1];
  int intersectx = - max( delta.x, deltaShape.x ) + min( delta.x + size.x, deltaShape.x + sizeShape.x );
  int intersecty = - max( delta.y, deltaShape.y ) + min( delta.y + size.y, deltaShape.y + sizeShape.y );
  if ( ( intersectx > 0 ) || ( intersecty > 0 ) ){
    T point;
    int direction = 0; // its southeast at first
    if ( delta.x > deltaShape.x ){
      point.x = deltaShape.x + sizeShape.x;
      direction += 1; // direction is changed to north
    }
    else
      point.x = deltaShape.x;
    if ( delta.y > deltaShape.y ){
      point.y = deltaShape.y + sizeShape.y;
      direction += 2; // direction is changed to west
    }
    else
      point.y = deltaShape.y;
    
    v.intersectionPoint = point;
    v.direction = direction;
    
  }
  if ( ( intersectx > 0 ) || ( intersecty > 0 ) ) return true;
  return false;
}

/**
 * \fn bool compareVertex( RefVertex< T > va, RefVertex< T > vb )
 * 
 * \brief This method compares two references vertex to be used in
 * the sort algorithm
 *
 * \param va, vb - References vertex that contains intersectionPoint 
 * and direction associated
 *
 * \return True if the x-axis of va is smaller than x-axis of vb and 
 * false otherwise.
 */
template <typename T>
bool
Blocks< T >::compareVertex(RefVertex< T > va, RefVertex< T > vb){
  return va.intersectionPoint.y < vb.intersectionPoint.y;
}


/** @} */ //end of group class.
