/**
 * \file polygons.hpp
 *
 * \brief This file contains the prototype of shape polygons.
 *
 * \author 
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date September 2019
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

#ifndef POLYGONS_HPP
#define POLYGONS_HPP

#include <iostream> //std::cout, std::endl
//#include "shape.hpp" //Shape
#include <vector> //std::vector
#include <math.h> //cos, sin
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

#define PI 3.1415 ///< Defining value to PI

using namespace std;
using namespace cv;

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Polygons
 *
 * \brief This class implements the Shape TAD to represent structure
 * of fovea.
 *
 * \tparam T - Generic representation for type cv::Point
 */
template < typename T > // cv::Point 
class Polygons : public Shape< T > {
public:
  //
  // Methods
  //  
  /**
   * \fn Polygons( vector< T > boundingBox, int nVertices )
   *
   * \brief Constructor default.
   * This method initialize all vertices of the shape.
   *
   * \param boundingBox - Vector containing 2 positions with tuple
   * information the limits of rectangular region ( delta and size )
   * to create the fovea
   *        nVertices - Number of vertices
   */
  Polygons( vector< T > boundingBox, int nVertices );
  
  /**
   * \fn ~Polygons()
   *
   * \brief Destructor default
   */
  ~Polygons();
  
  /**
   * \fn void printVertices()
   *
   * \brief Pure virtual method caracterize this class like abstract.
   * This method print all vertices of the shape.
   */
  void printVertices();
    
};

#endif

/**
 * \fn Polygons( vector< T > boundingBox, int nVertices )
 *
 * \brief Constructor default.
 * This method initialize all vertices of the shape.
 *
 * \param boundingBox - Vector containing 2 positions with tuple
 * information the limits of rectangular region ( delta and size )
 * to create the fovea
 *        nVertices - Number of vertices
 */
template< typename T >
Polygons< T >:: Polygons( vector< T > boundingBox, int nVertices ) : Shape< T >( boundingBox ){
  T delta = boundingBox[0];
  T size = boundingBox[1];
  vector< T > v;
  for ( int i = 0; i < nVertices; i++ ){
    float angle = 2*PI*i/nVertices;
    //v.push_back( T((int)cos(angle), (int)sin(angle)) );
    // Coordenada x = ( Sx * (m + 1) + 2Dx )/2, onde m = cos(angle)
    // Coordenada y = ( Sy * (n + 3) + 2Dy )/2, onde n = sin(angle)
    // Using conversion between dimensions -1 to 1 (model) to delta and size (image)
    int m = (int)cos(angle); int n = (int)sin(angle);
    float x_axis = (size.x + (m + 1) + 2*( delta.x ))/2;
    float y_axis = (size.y + (n + 3) + 2*( delta.y ))/2;
    v.push_back( T( x_axis, y_axis ) );
  }
  this->updateVertices( v );
}

/**
 * \fn ~Polygons()
 *
 * \brief Destructor default
 */
template< typename T >
Polygons< T >::~Polygons(){
  vector< T >().swap( this->vertices ); // Free the memory
}

/**
 * \fn void printVertices()
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * This method print all vertices of the shape.
 */
template< typename T >
void 
Polygons< T >::printVertices(){
  for ( int i = 0; i < this->vertices.size(); i++ )
    cout << "(x, y) = (" << this->vertices[i].x << ", " << this->vertices[i].y << ")" << endl;
}

/** @} */ //end of group class.
