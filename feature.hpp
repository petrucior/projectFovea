/**
 * \file feature.hpp
 *
 * \brief This file contains the prototype of feature TAD, responsable to
 * define extraction and computation features in the fovea.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
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

#ifndef FEATURE_HPP
#define FEATURE_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Feature
 *
 * \brief This class implements the Feature TAD to represent feature
 * of fovea with generic types ( feature and levels ).
 *
 * \tparam T - Generic representation for type cv::Point
 * \tparam K - Generic representation for type cv::ORB, cv::SURF
 */
template < typename T, typename K > // cv::Point / cv::Ptr<ORB>, cv::Ptr<SURF>
class Feature {  
public:
  //
  // Methods
  //
  /**
   * \fn Feature()
   *
   * \brief Constructor default.
   *
   */
  Feature();
  
  /**
   * \fn ~Feature()
   *
   * \brief Destructor default
   */
  ~Feature();
  
  /**
   * \fn virtual void show() = 0;
   *
   * \brief Pure virtual method caracterize this class like abstract.
   */
  virtual void show() = 0;
    
};

#endif

/**
 * \fn Feature()
 *
 * \brief Constructor default.
 */
template <typename T, typename K>
Feature< T, K >::Feature(){ 
  // No need to implement this method
}

/**
 * \fn ~FoveatedLevel()
 *
 * \brief Destructor default
 */
template <typename T, typename K>
Feature< T, K >::~Feature(){
  // No need to implement this method!
}

/** @} */ //end of group class.
