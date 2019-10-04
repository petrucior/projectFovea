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
//#include "level.hpp" //std::vector< Level >
//#include "feature.hpp"
#include "fovea.hpp"
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

  /**
   * \fn cv::Mat multifoveatedImage( cv::Mat img )
   *
   * \brief This function builds an image with multiples focus.
   *
   * \param img - Image to be multifoveated
   *
   * \return Image multifoveated created by multiples focus
   */
  cv::Mat multifoveatedImage( cv::Mat img );
  
private:
  //
  // Methods
  //
  
  //
  // Attributes
  //
  std::vector< Fovea< T >* > foveas; ///< Pointers to foveas
  // Parameters
  int m; ///< Number levels of fovea
  T w; ///< Size of levels
  T u; ///< Size of image
  T f; ///< Position (x, y) to build the fovea
  
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
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
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
  // Keeping values
  this->m = m;
  this->w = w;
  this->u = u;
  this->f = f;
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
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
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
  // Keeping values
  this->m = m;
  this->w = w;
  this->u = u;
  this->f = f;
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

/**
 * \fn cv::Mat multifoveatedImage( cv::Mat img )
 *
 * \brief This function builds an image with multiples focus.
 *
 * \param img - Image to be multifoveated
 *
 * \return Image multifoveated created by multiples focus
 */
template <typename T>
cv::Mat 
Multifovea< T >::multifoveatedImage( cv::Mat img ){
  std::vector<cv::Scalar> colors;
  srand (time(NULL)); // Initialize random seed
  for (int i = 0; i < foveas.size(); i++){
    int r = rand() % 256; // 0 - 255
    int g = rand() % 256; // 0 - 255
    int b = rand() % 256; // 0 - 255
    colors.push_back(cv::Scalar(b, g, r));
  }
  cv::Mat imgMultifoveated = img.clone();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m) // Schedule(static, m) keeps the order
#endif
  for ( int k =  0; k < this->m; k++ ){ // Levels
    for ( int focus = 0; focus < this->foveas.size(); focus++ ){ // foveas
      Level< T > level = (this->foveas[focus])->getLevelFromFovea( k );
      cv::Mat imgLevel = level.getLevel( img );
      // Mapping levels to foveated image
      std::vector< T > mapLevel2Image = (foveas[focus])->getMapLevel2ImageFoveaByLevel( k );
      T initial = mapLevel2Image[0];
      T final = mapLevel2Image[1];
#ifdef DEBUG
      std::cout << "(xi, yi) = (" << initial.x << ", " << initial.y << ")" << std::endl;
      std::cout << "(xf, yf) = (" << final.x << ", " << final.y << ")" << std::endl;
#endif
      cv::Rect roi = cv::Rect( initial.x, initial.y, final.x - initial.x, final.y - initial.y );
      if ( k < m ){ // Copying levels to foveated image
	resize( imgLevel, imgLevel, cv::Size(final.x - initial.x, final.y - initial.y), 0, 0, CV_INTER_LINEAR );
	imgLevel.copyTo( imgMultifoveated( roi ) );
      }
      else
	imgLevel.copyTo( imgMultifoveated( roi ) );
      
      //
      // Text and rectangles
      //
      /*
      char buffer[50];
      sprintf(buffer, "Fovea %d", focus);
      putText(imgMultifoveated, buffer, T(10, (10*focus)+10), cv::FONT_HERSHEY_SIMPLEX, 0.25, colors[focus], 1, 1);
      // Paint rectangle in each level
      cv::rectangle(imgMultifoveated, cv::Point(initial.x, initial.y), cv::Point(final.x - 1, final.y - 1), colors[focus]);
      */
    }
  }
  
  return imgMultifoveated;
}

/** @} */ //end of group class.