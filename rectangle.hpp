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
#include "shape.hpp" //Shape
#include <vector> //std::vector
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

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
 * \tparam T - Generic representation for type cv::Point
 */
template < typename T > // cv::Point 
class Rectangle : public Shape {
public:
  //
  // Methods
  //
  /**
   * \fn Rectangle()
   *
   * \brief Constructor default
   */
  Rectangle();
  
  /**
   * \fn void defVertices( std::vector< T > boundingBox )
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method initialize all vertices of the shape.
   *
   * \param boundingBox - Vector containing 2 positions with tuple
   * information the limits of rectangular region ( delta and size )
   * to create the fovea.
   */
  void defVertices( std::vector< T > boundingbox );
  
  //
  // It's necessary to make an overwrite in the function testIntersection
  // because this class implement a classic approach of fovea model
  //

  /**
   * \fn bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection )
   *
   * \brief This method checks if the shape intersect the region delimited
   * by vertices and return the point of intersection.
   *
   * \param shape - Shape analyzed
   * \param pointIntersection - List of positions (x, y) that intersect the shape
   *
   * \return True if shape intersect and false otherwise.
   */
  //bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection );
  
private:
  //
  // Attributes
  //
  T delta;
  T size;

};

#endif

/**
 * \fn Rectangle()
 *
 * \brief Constructor default
 */
Rectangle:: Rectangle(){
  // No need to implement!
}
  

/**
 * \fn void defVertices( std::vector< T > boundingBox )
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * This method initialize all vertices of the shape.
 *
 * \param boundingBox - Vector containing 2 positions with tuple
 * information the limits of rectangular region ( delta and size )
 * to create the fovea.
 */
void defVertices( std::vector< T > boundingbox ){
  this->delta = boundingBox[0];
  this->size = boundingBox[1];
}

/**
 * \fn bool intersectionShape( Shape& shape, std::vector< T >& pointIntersection )
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
Rectangle::intersectionShape( Shape& shape, std::vector< T >& pointIntersection ){
  // To do
}*/

/** @} */ //end of group class.
