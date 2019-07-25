/**
 * \file fovea.hpp
 *
 * \brief This file contains the prototype of strategies:
 * (1) extracting features of levels and (2) extracting 
 * octaves of features on levels.
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
 * This file is part of projectFoveaCuda software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 2 of the License, or (at your option) any later
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FOVEA_HPP
#define FOVEA_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#include "level.hpp" //std::vector< Level >
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Fovea
 *
 * \brief This class implements the Fovea TAD with generic type
 *
 * \tparam T - Generic representation for type cv::Point
 */
template< typename T > // cv::Point
class Fovea{
public:
  //
  // Methods
  //
  /**
   * \fn Fovea(int m, T w, T u, T f)
   *
   * \brief Constructor default of fovea class. THis method was 
   * built 
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   */
  Fovea(int m, T w, T u, T f);

private:
  //
  // Methods
  //
  /**
   * \fn void checkParameters( int m, T w, T u, T f )
   *
   * \brief This method check the parameters to build a fovea
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   */
  void checkParameters( int m, T w, T u, T f );
  
  //
  // Attributes
  //
  std::vector< Level > levels; ///< List of levels
  int m; ///< Number levels of fovea
  T w; ///< Size of levels
  T u; ///< Size of image
  T f; ///< Position (x, y) to build the fovea
  
};

#endif


/**
 * \fn Fovea(int m, T w, T u, T f);
 *
 * \brief Constructor default of fovea class
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 */
Fovea::Fovea(int m, T w, T u, T f){
  this->checkParameters( m, w, u, f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m; k++ ){
    Level l = new Level( k, m, w, u, f );
    levels.push_back( l );
  }
}

/**
 * \fn void checkParameters( int m, T w, T u, T f )
 *
 * \brief This method check the parameters to build a fovea
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 */
void 
Fovea::checkParameters( int m, T w, T u, T f ){
  // Verify if there is the minimun one level
  assert( m >= 1 );
  // Verify if w is bigger than zero and lower than image size
  assert( ( w.x > 0 ) && ( w.x < u.x ) );
  assert( ( w.y > 0 ) && ( w.y < u.y ) );
  // Verify if u is bigger than zero
  assert( ( u.x > 0 ) && ( u.y > 0 ) );
}


/** @} */ //end of group class.
