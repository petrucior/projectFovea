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
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "level.hpp" //std::vector< Level >
#include "feature.hpp"
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

#define MRMF 0 ///< Identify the use of MRMF method
#define MMF 1 ///< Identify the use of MMF method

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
   * \fn Fovea(cv::Mat img, int m, T w, T f)
   *
   * \brief Constructor default of fovea class.
   * This constructor create a fovea associeted to an image.
   * This approcah can be applied to associate the fovea structure
   * to unique image. 
   *
   * \param img - Image to be foveated
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param f - Position (x, y) of the fovea
   */
  Fovea(cv::Mat img, int m, T w, T f);
  
  /**
   * \fn Fovea(int m, T w, T u, T f)
   *
   * \brief Constructor default of fovea class.
   * This constructor create a fovea associated to image parameters.
   * This approach can be used and explored by robots or systems
   * which doesn't depends of the different parameters or device.  
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) of the fovea
   */
  Fovea(int m, T w, T u, T f);
  
  /**
   * \fn ~Fovea()
   *
   * \brief Destructor default of fovea class.
   */
  ~Fovea();
  
  /**
   * \fn void updateFovea(int m, T w, T u, T f)
   *
   * \brief This method update the fovea structure
   * and can be thought in update only what need to 
   * be changed, in other words, the update occors
   * only in shape of level. 
   * It is necessary to be implemented!
   *
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) of the fovea
   */
  void updateFovea(int m, T w, T u, T f);

  /**
   * \fn bool foveatedFeatures( cv::Mat img, Feature feature, int code )
   * \fn bool computeAndExtractFeatures( cv::Mat img, Feature< T > feature, int code )
   *
   * \brief This method compute and extract features 
   * of foveated structure using MRMF or MMF.
   *
   * \param img - Image to be foveated
   * \param code - This code indicates which method
   * to foveation will be used. If code is zero, then
   * MRMF is chosen, otherwise MMF.
   *
   * \return True if was done computed and extracted
   * features and False otherwise.
   */
  bool foveatedFeatures( cv::Mat img, int feature, int code );
  
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
  std::vector< Level< T > > levels; ///< List of levels
  int m; ///< Number levels of fovea
  T w; ///< Size of levels
  T u; ///< Size of image
  T f; ///< Position (x, y) to build the fovea
};

#endif

/**
 * \fn Fovea(cv::Mat img, int m, T w, T f)
 *
 * \brief Constructor default of fovea class.
 *
 * \param img - Image to be foveated
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param f - Position (x, y) of the fovea
 */
template <typename T>
Fovea< T >::Fovea(cv::Mat img, int m, T w, T f){
  T u = T( img.cols, img.rows );
  this->checkParameters( m, w, u, f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m; k++ ){
    Level< T > l( k, m, w, u, f );
    levels.push_back( l );
  }
}

/**
 * \fn Fovea(int m, T w, T u, T f)
 *
 * \brief Constructor default of fovea class.
 * This constructor create a fovea associated to image parameters.
 * This approach can be used and explored by robots or systems
 * which doesn't depends of the different parameters or device.  
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) of the fovea
 */
template <typename T>
Fovea< T >::Fovea(int m, T w, T u, T f){
  this->checkParameters( m, w, u, f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m; k++ ){
    Level< T > l( k, m, w, u, f );
    levels.push_back( l );
  }
} 

/**
 * \fn ~Fovea()
 *
 * \brief Destructor default of fovea class.
 */
template <typename T>
Fovea< T >::~Fovea(){
  std::vector< Level< T > >().swap( this->levels ); // Free the memory
}

/**
 * \fn void updateFovea(int m, T w, T u, T f)
 *
 * \brief This method update the fovea structure
 * and can be thought in update only what need to 
 * be changed, in other words, the update occors
 * only in shape of level. 
 * It is necessary to be implemented!
 *
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) of the fovea
 */
template <typename T>
void 
Fovea< T >::updateFovea(int m, T w, T u, T f){
  this->checkParameters( m, w, u, f );
  if ( levels.size() > 0 ){
    // cleaning vector levels, but conserving the first level
    levels.erase( levels.begin() + 1, levels.end() );
    levels.shrink_to_fit();
  }
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = levels.size(); k < m; k++ ){
    Level< T > l( k, m, w, u, f );
    levels.push_back( l );
  }
}

/**
 * \fn bool foveatedFeatures( cv::Mat img, Feature feature, int code )
 *
 * \brief This method compute and extract features 
 * of foveated structure using MRMF or MMF.
 *
 * \param img - Image to be foveated
 * \param code - This code indicates which method
 * to foveation will be used. If code is zero, then
 * MRMF is chosen, otherwise MMF.
 *
 * \return True if was done computed and extracted
 * features and False otherwise.
 */
template <typename T>
bool 
Fovea< T >::foveatedFeatures( cv::Mat img, int feature, int code ){
  if ( code == MRMF ){
    std::cout << "MRMF actived" << std::endl;
    //int feature = _KAZE_;
    Feature< cv::Point, int > activedMRMF( img, levels, feature );
    activedMRMF.show( img, levels );
  }
  if ( code == MMF ){
    std::cout << "MMF actived" << std::endl;
  }
  return true;
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
template <typename T>
void
Fovea< T >::checkParameters( int m, T w, T u, T f ){
  // Verify if there is the minimun one level
  assert( m >= 1 );
  // Verify if w is bigger than zero and lower than image size
  assert( ( w.x > 0 ) && ( w.x < u.x ) );
  assert( ( w.y > 0 ) && ( w.y < u.y ) );
  // Verify if u is bigger than zero
  assert( ( u.x > 0 ) && ( u.y > 0 ) );
}


/** @} */ //end of group class.
