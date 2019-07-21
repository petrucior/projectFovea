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
   * \fn Fovea( )
   *
   * \brief Constructor default of fovea class
   */
  Fovea();

private:
  //
  // Attributes
  //
  int m; ///< Number levels of fovea
  T w; ///< Size of levels
  T u; ///< Size of image
  T f; ///< Position (x, y) to build the fovea
  
}

#endif

/** @} */ //end of group class.
