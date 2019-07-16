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


/**
 * \class Shape
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea with generic type.
 */
Template < typename T > // cv::Point
class Shape {
public:
  //
  // Methods
  //
  /**
   * \fn virtual void defVertices( int m, T w, T u, T f ) = 0
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method initialize all vertices of the shape.
   *
   * \param m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) to build the fovea
   */
  virtual void defVertices( int m, T w, T u, T f ) = 0;
  
  /**
   * \fn std::vector< T > boundingBox( int k, int m, T w, T u, T f )
   *
   * \brief This method return the bounding box delimiting
   * the region where will be created the shape.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) to build the fovea
   *
   * \return Vector containing 2 positions with tuple 
   * information the limits of rectangular region ( delta and size ).
   */
  std::vector< T > boundingBox( int k, int m, T w, T u, T f );
  
  /**
   * \fn virtual bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection )
   *
   * \brief This method checks if the shape intersect the region delimited
   * by vertices and return the point of intersection.
   *
   * \param shape - Shape analyzed
   *        pointIntersection - List of positions (x, y) that intersect the shape
   *
   * \return True if shape intersect and false otherwise.
   */
  virtual bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection );    
  
private:
  //
  // Methods
  //
  /**
   * \fn T getDelta( int k, int m, T w, T u, T f )
   *
   * \brief Calculates the initial pixel to build MMF.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) to build the fovea
   *
   * \return Return the initial pixel on the both axis of level k to build MMF.
   */
  T getDelta( int k, int m, T w, T u, T f );

  /**
   * \fn T getSize( int k, int m, T w, T u )
   *
   * \brief Calculates the final pixel to build MMF.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *
   * \return Return the final pixel on the both axis of level k to build MMF.
   */
  T getSize( int k, int m, T w, T u );

  /**
   * \fn bool distance( T vertexA, T vertexB, T point )
   *
   * \brief Calculates an equation of the line and computes the value of distance
   * between the point and the line.
   *
   * \param vertexA - First vertex to build the line
   *        vertexB - Second vertex to build the line
   *        point - Point to be analyzed
   *
   * \return Return true if distance is bigger or equal to zero, in other words,
   * the point is up or on the right of the line and false otherwise ( left side ) 
   */
  bool distance( T vertexA, T vertexB, T point );
  
  
  //
  // Attributes
  //
  std::vector< T > vertices; // vertices of shape
  int _m; // Number levels of fovea
  T _w; // Size of levels
  T _u; // Size of image
  T _f; // Position (x, y) to build the fovea

};

#endif

/**
 * \fn std::vector< T > boundingBox( int k, int m, T w, T u, T f )
 *
 * \brief This method return the bounding box delimiting
 * the region where will be created the shape.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) to build the fovea
 *
 * \return Vector containing 2 positions with tuple 
 * information the limits of rectangular region ( delta and size ).
 */
std::vector< T > 
Shape::boundingBox( int k, int m, T w, T u, T f ){
  T delta = getDelta( k, m, w, u, f );
  T size = getSize( k, m, w, u );
  std::vector< T > _boundingBox;
  _boundingBox.push_back( delta );
  _boundingBox.push_back( size );
  return _boundingBox;
}

/**
 * \fn virtual bool intersectionShape( Shape& shape, T pointIntersection )
 *
 * \brief This method checks if the shape intersect the region delimited
 * by vertices and return the point of intersection.
 *
 * \param shape - Shape analyzed
 *        pointIntersection - Position (x, y) that intersect the shape
 *
 * \return True if shape intersect and false otherwise.
 */
virtual bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection ){
  T vertexShape, vertexA, vertexB;
  bool insideStruct;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int i = 0; i < shape->vertices.size(); i++ ){ // Shape vertices iterator
    vertexShape = shape->vertices[i];
    insideStruct = true;
    for ( int j = 0; (( j < vertices.size() ) && ( insideStruct != false )); j++ ){ // Current shape vertices iterator
      vertexA = vertices[j];
      vertexB = vertices[0];
      // Only the last vertex is not updated
      if ( j < vertices.size() - 1 )
	vertexB = vertices[j+1];
      
      if ( distance( vertexA, vertexB, vertexShape ) == false )
	insideStruct &= false;
      
    }
    // Insert point of intersection between shapes
    if ( insideStruct == true )
      pointIntersection.push_back( vertexShape );
  }
    
}


/**
 * \fn T getDelta( int k, int m, T w, T u, T f )
 *
 * \brief Calculates the initial pixel to build MMF.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) to build the fovea
 *
 * \return Return the initial pixel on the both axis of level k to build MMF.
 */
T 
Shape::getDelta( int k, int m, T w, T u, T f ){
  int dx = int( k * ( u.x - w.x + ( 2 * f.x ) ) )/ ( 2 * m );
  int dy = int( k * ( u.y - w.y + ( 2 * f.y ) ) )/ ( 2 * m );
#ifdef DEBUG
  std::cout << "Delta: ( " << dx << ", " << dy << " ) " << std::endl;  
#endif
  return T( dx, dy );
}

/**
 * \fn T getSize( int k, int m, T w, T u )
 *
 * \brief Calculates the final pixel to build MMF.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *
 * \return Return the final pixel on the both axis of level k to build MMF.
 */
T 
Shape::getSize( int k, int m, T w, T u ){
  int sx = ((m * u.x) + (w.x * k) - (k * u.x)) / m;
  int sy = ((m * u.y) + (w.y * k) - (k * u.y)) / m;
#ifdef DEBUG
  std::cout << "Size: ( " << sx << ", " << sy << " ) " << std::endl;  
#endif
  return T( sx, sy );
}

/**
 * \fn bool distance( T vertexA, T vertexB, T point )
 *
 * \brief Calculates an equation of the line and computes the value of distance
 * between the point and the line.
 *
 * \param vertexA - First vertex to build the line
 *        vertexB - Second vertex to build the line
 *        point - Point to be analyzed
 *
 * \return Return true if distance is bigger or equal to zero, in other words,
 * the point is up or on the right of the line and false otherwise ( left side ) 
 */
bool 
Shape::distance( T vertexA, T vertexB, T point ){
  float a = vertexA.y - vertexB.y;
  float b = vertexB.x - vertexA.x;
  float c = (vertexA.x * vertexB.y) - (vertexB.x * vertexA.y);
  float dist = ( a * point.x + b * point.y + c ) / ( sqrt ( pow( a, 2.0 ) + pow( b, 2.0 ) ) ); 
  if ( dist >= 0.0 ) return true;
  return false;
}

