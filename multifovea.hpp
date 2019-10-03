/**
 * \file multifovea.hpp
 *
 * \brief This file contains the prototype of strategies to 
 * multifoveation: (1) reexecution, (2) pixel by pixel approach, 
 * (3) bitmap and (4) blocks sending.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
 *
 * \version 0.1
 * \date October 2019
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

#ifndef MULTIFOVEA_HPP
#define MULTIFOVEA_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "fovea.hpp"
#include "level.hpp" //std::vector< Level >
#include "feature.hpp"
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

// Setting Multifovea Approaches
#define REEXECUTION 0 ///< Identify the use of reexecution approach
#define PIXELBYPIXEL 1 ///< Identify the use of pixel-by-pixel processing approach
#define BITMAP 2 ///< Identify the use of bitmap approach
#define BLOCKS 3 ///< Identify the use of block-based processing approach

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class MultiFovea
 *
 * \brief This class implements the Multifovea TAD with generic type
 *
 * \tparam T - Generic representation for type cv::Point
 *
 * \note
 * This class works only for similar structures, in other words, use 
 * the same parameters for multiple fovea.
 */
template< typename T > // cv::Point
class Multifovea{
public:
  //
  // Methods
  //
  /**
   * \fn Multifovea(cv::Mat img, int m, T w, std::vector< T > fs, int mode)
   *
   * \brief Constructor default of multifovea class.
   * This constructor create many foveas associeted to an image.
   * This approcah can be applied to associate the foveas structure
   * to unique image. 
   *
   * \param img - Image to be foveated
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param fs - Vector with positions (x, y) of the foveas
   * \param mode - identify the approach: reexecution, pixel-by-pixel, 
   * bitmap or block-based (see settings multifovea approaches )
   */
  Multifovea(cv::Mat img, int m, T w, std::vector< T > fs, int mode);
  
  /**
   * \fn Multifovea(int m, T w, T u, std::vector< T > fs, int mode)
   *
   * \brief Constructor default of multifovea class.
   * This constructor create many foveas associated to image parameters.
   * This approach can be used and explored by robots or systems
   * which doesn't depends of the different parameters or device.  
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param fs - Vector with positions (x, y) of the foveas
   * \param mode - identify the approach: reexecution, pixel-by-pixel, 
   * bitmap or block-based (see settings multifovea approaches )
   */
  Multifovea(int m, T w, T u, std::vector< T > fs, int mode);
  
  /**
   * \fn ~Multifovea()
   *
   * \brief Destructor default of multifovea class.
   */
  ~Multifovea();
  
private:
  //
  // Methods
  //
  
  //
  // Attributes
  //
  std::vector< Fovea< T >* > foveas; ///< Pointers to foveas

};

#endif

/**
 * \fn Multifovea(cv::Mat img, int m, T w, std::vector< T > fs, int mode)
 *
 * \brief Constructor default of multifovea class.
 * This constructor create many foveas associeted to an image.
 * This approcah can be applied to associate the foveas structure
 * to unique image. 
 *
 * \param img - Image to be foveated
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param fs - Vector with positions (x, y) of the foveas
 * \param mode - identify the approach: reexecution, pixel-by-pixel, 
 * bitmap or block-based (see settings multifovea approaches )
 */
template <typename T>
Multifovea< T >::Multifovea(cv::Mat img, int m, T w, std::vector< T > fs, int mode){
  Fovea< T > *fovea;
  switch ( mode ){
  case REEXECUTION:
    std::cout << "reexecution approach" << std::endl;
    for ( int f = 0; f < fs.size(); f++ ){
      fovea = new Fovea< T >( img, m, w, fs[f] );
      foveas.push_back( fovea );
    }
    break;
  case PIXELBYPIXEL:
    std::cout << "pixel by pixel processing approach" << std::endl;
    break;
  case BITMAP:
    std::cout << "bitmap approach" << std::endl;
    break;
  case BLOCKS:
    std::cout << "Block-based processing approach" << std::endl;
    break;
  default:
    std::cout << "There was not configured the mode" << std::endl;
    break;
  }

}


/**
 * \fn Multifovea(int m, T w, T u, std::vector< T > fs, int mode)
 *
 * \brief Constructor default of multifovea class.
 * This constructor create many foveas associated to image parameters.
 * This approach can be used and explored by robots or systems
 * which doesn't depends of the different parameters or device.  
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param fs - Vector with positions (x, y) of the foveas
 * \param mode - identify the approach: reexecution, pixel-by-pixel, 
 * bitmap or block-based (see settings multifovea approaches )
 */
template <typename T>
Multifovea< T >::Multifovea(int m, T w, T u, std::vector< T > fs, int mode){
  Fovea< T > *fovea;
  switch ( mode ){
  case REEXECUTION:
    std::cout << "reexecution approach" << std::endl;
    for ( int f = 0; f < fs.size(); f++ ){
      fovea = new Fovea< T >( m, w, u, fs[f] );
      foveas.push_back( fovea );
    }
    break;
  case PIXELBYPIXEL:
    std::cout << "pixel by pixel processing approach" << std::endl;
    break;
  case BITMAP:
    std::cout << "bitmap approach" << std::endl;
    break;
  case BLOCKS:
    std::cout << "Block-based processing approach" << std::endl;
    break;
  default:
    std::cout << "There was not configured the mode" << std::endl;
    break;
  }
  
}


/**
 * \fn ~Multifovea()
 *
 * \brief Destructor default of multifovea class.
 */
template <typename T>
Multifovea< T >::~Multifovea(){
  std::vector< Fovea< T > >().swap( this->foveas ); // Free the memory
}

/** @} */ //end of group class.
