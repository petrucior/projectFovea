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
#include <vector> //std::vector

/**
 * \class Rectangle
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea.
 */
Template < typename T > // cv::Point 
class Rectangle : public Shape {
public:
  //
  // Methods
  //
  /**
   * \fn void defVertices( int m, T w, T u, T f )
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method initialize all vertices of the shape.
   * \param m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) to build the fovea
   */
  void defVertices( int m, T w, T u, T f );

  //
  // It's necessary to make an overwrite in the function testIntersection
  // because this class implement a classic approach of fovea model
  //
  
};

#endif

/**
 * \fn void defVertices( int m, T w, T u, T f )
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * This method initialize all vertices of the shape.
 * \param m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) to build the fovea
 */
void 
Rectangle::defVertices( int m, T w, T u, T f ){
  vertices = boundingBox( int k, int m, T w, T u, T f );
}
