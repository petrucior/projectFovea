/**
 * \file shape.hpp
 *
 * \brief This file contains the prototype of shapes to build the 
 * fovea structure.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
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

#ifndef SHAPE_HPP
#define SHAPE_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#include <math.h> //pow, sqrt
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

using namespace std;

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Shape
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea with generic type.
 *
 * \tparam T - Generic representation for type cv::Point
 */
template < typename T > // cv::Point
class Shape {  
protected:  
  //
  // Attributes
  //
  vector< T > vertices; ///< vertices of shape ( clockwise )

public:
  //
  // Methods
  //
  /**
   * \fn Shape( vector< T > boundingBox )
   *
   * \brief Constructor default.
   * This method initialize all vertices of the shape.
   *
   * \param boundingBox - Vector containing 2 positions with tuple
   * information the limits of rectangular region ( delta and size ) 
   * to create the fovea
   * Format: [ (delta.x, delta.y), (size.x, delta.y), (size.x, size.y), (delta.x, size.y)]
   */
  Shape( vector< T > boundingBox ); //{ this->vertices = boundingBox; }

  /**
   * \fn ~Shape()
   *
   * \brief Destructor default
   */
  ~Shape(); //{ vector< T >().swap( this->vertices ); } // Free the memory

  /**
   * \fn void updateVertices( vector< T > shapeVertices )
   *
   * \brief Update vertices of shape class.
   *
   * \param shapeVertices - all vertices to be updated.
   */
  void updateVertices( vector< T > shapeVertices );

  /**
   * \fn virtual void printVertices() = 0
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method print all vertices of the shape.
   */
  virtual void printVertices() = 0;
  
  /**
   * \fn virtual vector< T > getVertices()
   *
   * \brief This function informs the vertices of shape
   *
   * \return Return the shape constructed by level
   */
  virtual vector< T > getVertices();

  /**
   * \fn virtual bool intersectionShape( T point )
   *
   * \brief This method checks if the shape is intersected by point.
   *
   * \param point - Point (x, y) analyzed
   *
   * \return True if shape is intersected and false otherwise.
   */
  virtual bool intersectionShape( T  point );
  

  /**
   * \fn virtual bool intersectionShape( const Shape< T >& shape, vector< T >& pointIntersection )
   *
   * \brief This method checks if the shape intersect the region delimited
   * by vertices and return the point of intersection.
   *
   * \param shape - Shape analyzed
   * \param pointIntersection - List of positions (x, y) that intersect the shape
   *
   * \return True if shape intersect and false otherwise.
   */
  virtual bool intersectionShape( const Shape< T >& shape, vector< T >& pointIntersection );    
  
private:
  //
  // Methods
  //
  /**
   * \fn bool distance( T vertexA, T vertexB, T point )
   *
   * \brief Calculates an equation of the line and computes the value of distance
   * between the point and the line.
   *
   * \param vertexA - First vertex to build the line
   * \param vertexB - Second vertex to build the line
   * \param point - Point to be analyzed
   *
   * \return Return true if distance is bigger or equal to zero, in other words,
   * the point is up or on the right of the line and false otherwise ( left side ) 
   */
  bool distance( T vertexA, T vertexB, T point );
  
};

#endif

/**
 * \fn Shape( vector< T > boundingBox )
 *
 * \brief Constructor default.
 * This method initialize all vertices of the shape (clockwise).
 *
 * \param boundingBox - Vector containing 2 positions with tuple
 * information the limits of rectangular region ( delta and size ) 
 * to create the fovea
 * Format: [ (delta.x, delta.y), (size.x, delta.y), (size.x, size.y), (delta.x, size.y)]
 */
template <typename T>
Shape< T >::Shape( vector< T > boundingBox ){
  this->vertices = boundingBox;
}

/**
 * \fn ~Shape()
 *
 * \brief Destructor default
 */
template <typename T>
Shape< T >::~Shape(){
  vector< T >().swap( this->vertices ); // Free the memory
}

/**
 * \fn void updateVertices( vector< T > shapeVertices )
 *
 * \brief Update vertices of shape class.
 *
 * \param shapeVertices - all vertices to be updated.
 */
template <typename T>
void 
Shape< T >::updateVertices( vector< T > shapeVertices ){
  this->vertices = shapeVertices;
  //for ( int i = 0; i < shapeVertices.size(); i++ )
  //  cout <<  shapeVertices[i].x << ", " << shapeVertices[i].y << endl;
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
Shape< T >::getVertices(){
  return this->vertices;
}

/**
 * \fn virtual bool intersectionShape( T point )
 *
 * \brief This method checks if the shape is intersected by point.
 *
 * \param point - Point (x, y) analyzed
 *
 * \return True if shape is intersected and false otherwise.
 */
template <typename T>
bool
Shape< T >::intersectionShape( T  point ){
  T vertexA, vertexB;
  bool insideStruct = true;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int i = 0; (( i < this->vertices.size() ) && ( insideStruct != false )); i++ ){ // Shape vertices iterator
    vertexA = this->vertices[i];
    vertexB = this->vertices[0];
    // Only the last vertex is not updated
    if ( i < this->vertices.size() - 1 )
      vertexB = this->vertices[i+1];
    
    if ( distance( vertexA, vertexB, point ) == false )
      insideStruct = false;
  }
  return insideStruct;
  
}

/**
 * \fn virtual bool intersectionShape( const Shape< T >& shape, vector< T >& pointIntersection )
 *
 * \brief This method checks if the shape intersect the region delimited
 * by vertices and return the point of intersection.
 *
 * \param shape - Shape analyzed
 * \param pointIntersection - Position (x, y) that intersect the shape
 *
 * \return True if shape intersect and false otherwise.
 */
template < typename T >
bool 
Shape< T >::intersectionShape( const Shape< T >& shape, vector< T >& pointIntersection ){
  T vertexShape, vertexA, vertexB;
  bool insideStruct;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int i = 0; i < shape.vertices.size(); i++ ){ // Shape vertices iterator
    vertexShape = shape.vertices[i];
    insideStruct = true;
    for ( int j = 0; (( j < this->vertices.size() ) && ( insideStruct != false )); j++ ){ // Current shape vertices iterator
      vertexA = this->vertices[j];
      vertexB = this->vertices[0];
      // Only the last vertex is not updated
      if ( j < this->vertices.size() - 1 )
	vertexB = this->vertices[j+1];
      
      if ( distance( vertexA, vertexB, vertexShape ) == false )
	insideStruct &= false;
      
    }
    // Insert point of intersection between shapes
    if ( insideStruct == true )
      pointIntersection.push_back( vertexShape );
    
  }
  if ( pointIntersection.size() > 0 ) return true;
  else return false;
}

/**
 * \fn bool distance( T vertexA, T vertexB, T point )
 *
 * \brief Calculates an equation of the line and computes the value of distance
 * between the point and the line.
 *
 * \param vertexA - First vertex to build the line ( clockwise )
 * \param vertexB - Second vertex to build the line ( clockwise )
 * \param point - Point to be analyzed
 *
 * \return Return true if distance is bigger or equal to zero, in other words,
 * the point is up or on the right of the line and false otherwise ( left side ) 
 */
template <typename T>
bool 
Shape< T >::distance( T vertexA, T vertexB, T point ){
  double a = vertexA.y - vertexB.y;
  double b = vertexB.x - vertexA.x;
  double c = (vertexA.x * vertexB.y) - (vertexB.x * vertexA.y);
  double dist = ( a * point.x + b * point.y + c ) / ( sqrt ( pow( a, 2.0 ) + pow( b, 2.0 ) ) ); 
  if ( dist >= 0.0 ) return true;
  return false;
}


/** @} */ //end of group class.
