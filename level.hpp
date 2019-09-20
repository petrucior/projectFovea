/**
 * \file level.hpp
 *
 * \brief This file contains the prototype of levels.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
 *
 * \version 0.1
 * \date July 2019
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

#ifndef LEVEL_HPP
#define LEVEL_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "shape.hpp" //Shape
#include "rectangle.hpp" // Shape< T >::Rectangle

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Level
 *
 * \brief This class implements the Level TAD to represent levels
 * of fovea with generic type.
 *
 * \tparam T - Generic representation for type cv::Point
 */
template < typename T > // cv::Point
class Level {
public:
  //
  // Methods
  //

  /**
   * \fn Level( int k, int m, T w, T u, T f )
   *
   * \brief Default constructor for Level class which will create 
   * a level of Rectangle shape.
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) to build the fovea
   */
  Level( int k, int m, T w, T u, T f );
  
  /**
   * \fn cv::Mat getLevel( cv::Mat img )
   *
   * \brief Method responsable to create an image with 
   * with level dimension
   *
   * \param img - Image will be foveated
   *
   * \return Image that represent the level
   */
  cv::Mat getLevel( cv::Mat img );
  
  /**
   * \fn std::vector< T > boundingBox( int k, int m, T w, T u, T f )
   *
   * \brief This method return the bounding box delimiting
   * the region where will be created the shape.
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) to build the fovea
   *
   * \return Vector containing 2 positions with tuple 
   * information the limits of rectangular region ( delta and size ).
   */
  std::vector< T > boundingBox( int k, int m, T w, T u, T f );
  
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
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) to build the fovea
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
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   *
   * \return Return the final pixel on the both axis of level k to build MMF.
   */
  T getSize( int k, int m, T w, T u );
  
  //
  // Attributes
  //
  Shape< T >* shape; ///< Shape of level
  int indexLevel; ///< Index of level
};

#endif

/**
 * \fn Level( int k, int m, T w, T u, T f )
 *
 * \brief Default constructor for Level class which will create 
 * a level of Rectangle shape.
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) to build the fovea
 */
template <typename T>
Level< T >::Level( int k, int m, T w, T u, T f ){
  indexLevel = k;
  std::vector< T > _boundingBox = this->boundingBox( k, m, w, u, f );
  Rectangle< T > *r = new Rectangle< T >( _boundingBox );
  shape = r;
}


/**
 * \fn cv::Mat getLevel( cv::Mat img )
 *
 * \brief Method responsable to create an image with 
 * with level dimension
 *
 * \param img - Image will be foveated
 *
 * \return Image that represent the level
 */
template <typename T>
cv::Mat
Level< T >::getLevel( cv::Mat img ){
  // It's necessary cut the image depending of vertices him
  std::vector< T > boundingBox = shape->getVertices();
  cv::Rect roi(boundingBox[0].x, boundingBox[0].y, boundingBox[1].x, boundingBox[1].y);
  cv::Mat croppedImage = img( roi );
  return croppedImage;
}

/**
 * \fn std::vector< T > boundingBox( int k, int m, T w, T u, T f )
 *
 * \brief This method return the bounding box delimiting
 * the region where will be created the shape.
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) to build the fovea
 *
 * \return Vector containing 2 positions with tuple 
 * information the limits of rectangular region ( delta and size ).
 */
template <typename T>
std::vector< T > 
Level< T >::boundingBox( int k, int m, T w, T u, T f ){
  T delta = getDelta( k, m, w, u, f ); ///< Delta is the upper left corner of bounding box
  T size = getSize( k, m, w, u ); ///< Size is the dimension between Delta and bottom right corner of bounding box
  std::vector< T > _boundingBox;  ///< Tuple vector containing delta and size
  _boundingBox.push_back( delta );
  _boundingBox.push_back( size );
  return _boundingBox;
}

/**
 * \fn T getDelta( int k, int m, T w, T u, T f )
 *
 * \brief Calculates the initial pixel to build MMF.
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) to build the fovea
 *
 * \return Return the initial pixel on the both axis of level k to build MMF.
 */
template <typename T>
T 
Level< T >::getDelta( int k, int m, T w, T u, T f ){
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
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 *
 * \return Return the final pixel on the both axis of level k to build MMF.
 */
template <typename T>
T 
Level< T >::getSize( int k, int m, T w, T u ){
  int sx = ((m * u.x) + (w.x * k) - (k * u.x)) / m;
  int sy = ((m * u.y) + (w.y * k) - (k * u.y)) / m;
#ifdef DEBUG
  std::cout << "Size: ( " << sx << ", " << sy << " ) " << std::endl;  
#endif
  return T( sx, sy );
}


/** @} */ //end of group class.
