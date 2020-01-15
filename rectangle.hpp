/**
 * \file rectangle.hpp
 *
 * \brief This file contains the prototype of shape rectangle of classic approach.
 *
 * \author 
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date June 2019
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

#ifndef RECTANGLE_HPP
#define RECTANGLE_HPP

#include <iostream> //std::cout, std::endl
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

/**
 * \class Rectangle
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea.
 *
 * \tparam T - Generic representation for type Point
 */
template < typename T > // cv::Point 
class Rectangle : public Shape< T > {
public:
  //
  // Methods
  //
  /**
   * \fn Rectangle( vector< T > boundingBox )
   *
   * \brief Constructor default.
   * This method initialize all vertices of the shape.
   *
   * \param boundingBox - Vector containing 2 positions with tuple
   * information the limits of rectangular region ( delta and size )
   * to create the fovea
   */
  Rectangle( vector< T > boundingBox );
  
  /**
   * \fn ~Rectangle()
   *
   * \brief Destructor default
   */
  ~Rectangle();
  
  /**
   * \fn virtual void printVertices()
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method print all vertices of the shape.
   */
  void printVertices();
  
  /**
   * \fn virtual vector< T > getVertices()
   *
   * \brief This function informs the vertices of shape
   *
   * \return Return the shape constructed by level
   */
  vector< T > getVertices();
  
  //
  // It's necessary to make an overwrite in the function testIntersection
  // because this class implement a classic approach of fovea model
  //
  
  /**
   * \fn bool intersectionShape( const Shape< T >& shape, const vector< T >& pointIntersection )
   *
   * \brief This method checks if the shape intersect the region delimited
   * by vertices and return the point of intersection.
   *
   * \param shape - Shape analyzed
   * \param pointIntersection - List of positions (x, y) that intersect the shape
   *
   * \return True if shape intersect and false otherwise.
   */
  //bool intersectionShape( const Shape< T >& shape, const vector< T >& pointIntersection );
  
private:
  //
  // Attributes
  //
  T delta;
  T size;

};

#endif

/**
 * \fn Rectangle( vector< T > boundingBox )
 *
 * \brief Constructor default.
 * This method initialize all vertices of the shape.
 *
 * \param boundingBox - Vector containing 2 positions with tuple
 * information the limits of rectangular region ( delta and size )
 * to create the fovea.
 * Format: [ (delta.x, delta.y), (size.x, delta.y), (size.x, size.y), (delta.x, size.y)]
 */
template <typename T>
Rectangle< T >::Rectangle( vector< T > boundingBox ) : Shape< T >( boundingBox ){
  delta = boundingBox[0];
  size = boundingBox[1];
  vector< T > v;
  v.push_back( T( delta.x, delta.y) );
  v.push_back( T( size.x, delta.y) );
  v.push_back( T( size.x, size.y) );
  v.push_back( T( delta.x, size.y) );
  this->updateVertices( v );
}

/**
 * \fn ~Rectangle()
 *
 * \brief Destructor default
 */
template <typename T>
Rectangle< T >::~Rectangle(){
  vector< T >().swap( this->vertices ); // Free the memory
}

/**
 * \fn void printVertices()
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * This method print all vertices of the shape.
 */
template <typename T>
void 
Rectangle< T >::printVertices(){
  cout << "delta: " << delta << ", size: " << size << endl;
}

/**
 * \fn virtual vector< T > getVertices()
 *
 * \brief This function informs the vertices of shape
 *
 * \return Return the shape constructed by level
 */
template <typename T>
vector< T >
Rectangle< T >::getVertices(){
  vector< T > boundingBox;
  boundingBox.push_back( delta );
  boundingBox.push_back( size );
  return boundingBox;
}

/**
 * \fn bool intersectionShape( const Shape< T >& shape, const vector< T >& pointIntersection )
 *
 * \brief This method checks if the shape intersect the region delimited
 * by vertices and return the point of intersection.
 *
 * \param shape - Shape analyzed
 * \param pointIntersection - List of positions (x, y) that intersect the shape
 *
 * \return True if shape intersect and false otherwise.
 */
/*bool 
Rectangle::intersectionShape( const Shape< T >& shape, const vector< T >& pointIntersection ){
  // To do
}*/


/** @} */ //end of group class.
