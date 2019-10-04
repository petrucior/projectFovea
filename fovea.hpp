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
   * \fn inline void fixFovea() 
   *
   * \brief Fix the fovea position: if fovea is outsite
   * image domain, snap it to the closest valid position 
   * independently for each coordinate
   */
  inline void fixFovea();
  
  /**
   * \fn void setFovea( T px )
   *
   * \brief Convert image points to axis fovea.
   *
   * \param px - Image points ( x, y ) 
   */
  void setFovea( T px );
  
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
  void updateFovea(T f);
  
  /**
   * \fn cv::Mat foveatedImage( cv::Mat img )
   *
   * \brief This function builds the focused image.
   *
   * \param img - Image to be foveated
   *
   * \return Image foveated created by levels
   */
  cv::Mat foveatedImage( cv::Mat img );
  
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

  /**
   * \fn Level< T > getLevelFromFovea( int k )
   *
   * \brief This method return levels from fovea
   *
   * \param k - Index for indicate which level will be returned
   *
   * \return Levels from fovea
   */
  Level< T > getLevelFromFovea( int k );

  /**
   * \fn std::vector< T > getMapLevel2Image( int k )
   *
   * \brief This method calculates the position of pixel on the 
   * level to image by level
   *
   * \param k - level of fovea
   *
   * \return Vector containing the boundingBox to map of level
   */
  std::vector< T > getMapLevel2Image( int k );
  
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

  /**
   * \fn T mapLevel2Image( int k, int m, T w, T u, T f, T px )
   *
   * \brief Calculates the position of pixel on the level to image.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        w - Size of levels
   *        u - Size of image
   *        f - Position (x, y) of the fovea
   *        px - Pixel (x, y) that we want to map.
   *
   * \return Return the position of pixel on the both axis to image.
   */
  T mapLevel2Image( int k, int m, T w, T u, T f, T px );
  
  //
  // Attributes
  //
  std::vector< Level< T > > levels; ///< List of levels
  Feature< T, int >* features = NULL; ///< Features
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
 * \fn inline void fixFovea() 
 *
 * \brief Fix the fovea position: if fovea is outsite
 * image domain, snap it to the closest valid position 
 * independently for each coordinate
 */
template <typename T>
inline void 
Fovea< T >::fixFovea(){
  f.x = MIN((u.x - w.x)/2, f.x);
  f.x = MAX((w.x - u.x)/2, f.x);
  f.y = MIN((u.y - w.y)/2, f.y);
  f.y = MAX((w.y - u.y)/2, f.y);
}  

/**
 * \fn inline void setFovea( T px )
 *
 * \brief Convert image points to axis fovea.
 *
 * \param px - Image points ( x, y ) 
 */
template <typename T>
void
Fovea< T >::setFovea( T px ){
  f.x = px.x - (u.x/2);
  f.y = px.y - (u.y/2);
  fixFovea();
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
Fovea< T >::updateFovea(T f){
  setFovea( f );
  this->checkParameters( m, w, u, this->f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m; k++ ){
    levels[k].updateLevel( m, w, u, this->f );
  }
}

/**
 * \fn cv::Mat foveatedImage( cv::Mat img )
 *
 * \brief This function builds the focused image.
 *
 * \param img - Image to be foveated
 *
 * \return Image foveated created by levels
 */
template <typename T>
cv::Mat 
Fovea< T >::foveatedImage( cv::Mat img ){
  cv::Mat imgFoveated = img.clone();
  std::vector< cv::KeyPoint > kp;
  std::vector< std::vector< cv::KeyPoint > > keypoints;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m+1) // Schedule(static, m+1) keeps the order
#endif
  for ( int k = 0; k < levels.size(); k++ ){ // Levels
    cv::Mat imgLevel = levels[k].getLevel( img );
    if ( features != NULL ){
      kp = features->getKeyPoints( k );
      for ( int i = 0; i < kp.size(); i++ ){
	cv::Point2f kpPos = mapLevel2Image( k, this->m, this->w, this->u, this->f, cv::Point2f( kp[i].pt.x, kp[i].pt.y ) );
	kp[i].pt.x = kpPos.x;
	kp[i].pt.y = kpPos.y;
      }
      keypoints.push_back( kp );
    }
    // Mapping levels to foveated image
    T initial = mapLevel2Image( k, this->m, this->w, this->u, this->f, T( 0, 0 ) ); 
    T final = mapLevel2Image( k, this->m, this->w, this->u, this->f, T( this->w.x, this->w.y ) );
#ifdef DEBUG
    std::cout << "(xi, yi) = (" << initial.x << ", " << initial.y << ")" << std::endl;
    std::cout << "(xf, yf) = (" << final.x << ", " << final.y << ")" << std::endl;
#endif
    cv::Rect roi = cv::Rect( initial.x, initial.y, final.x - initial.x, final.y - initial.y );
    if ( k < m ){ // Copying levels to foveated image
      resize( imgLevel, imgLevel, cv::Size(final.x - initial.x, final.y - initial.y), 0, 0, CV_INTER_LINEAR );
      imgLevel.copyTo( imgFoveated( roi ) );
    }
    else
      imgLevel.copyTo( imgFoveated( roi ) );
      
    // Paint rectangle in each level
    cv::rectangle(imgFoveated, cv::Point(initial.x, initial.y), cv::Point(final.x - 1, final.y - 1), cv::Scalar(255, 255, 255));

  }
  
  if ( keypoints.size() != 0 ){
    for ( int i = 0; i < keypoints.size(); i++ )
      cv::drawKeypoints( imgFoveated, keypoints[i], imgFoveated, cv::Scalar::all(-1), cv::DrawMatchesFlags::DEFAULT );
  }
  
  features = NULL;
  return imgFoveated;
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
    features = new Feature< T, int >( img, levels, feature );
#ifdef DEBUG
    features->show( img, levels );
#endif
  }
  if ( code == MMF ){
    std::cout << "MMF actived" << std::endl;
  }
  return true;
}

/**
 * \fn Level< T > getLevelFromFovea( int k )
 *
 * \brief This method return levels from fovea
 *
 * \param k - Index for indicate which level will be returned
 *
 * \return Levels from fovea
 */
template <typename T>
Level< T > 
Fovea< T >::getLevelFromFovea( int k ){
  return levels[k];
}

/**
 * \fn std::vector< T > getMapLevel2Image( int k )
 *
 * \brief This method calculates the position of pixel on the 
 * level to image by level
 *
 * \param k - level of fovea
 *
 * \return Vector containing the boundingBox to map of level
 */
template <typename T>
std::vector< T > 
Fovea< T >::getMapLevel2Image( int k ){
  std::vector< T > mapImage;
  mapImage.push_back( this->mapLevel2Image( k, this->m, this->w, this->u, this->f, T( 0, 0 ) ) );
  mapImage.push_back( this->mapLevel2Image( k, this->m, this->w, this->u, this->f, T( this->w.x, this->w.y ) ) );
  return mapImage;
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
  // Keeping values
  this->m = m;
  this->w = w;
  this->u = u;
  this->f = f;
}

/**
 * \fn T mapLevel2Image( int k, int m, T w, T u, T f, T px )
 *
 * \brief Calculates the position of pixel on the level to image.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        w - Size of levels
 *        u - Size of image
 *        f - Position (x, y) of the fovea
 *        px - Pixel (x, y) that we want to map.
 *
 * \return Return the position of pixel on the both axis to image.
 */
template <typename T>
T
Fovea< T >::mapLevel2Image( int k, int m, T w, T u, T f, T px ){
  int _px = ( (k * w.x) * (u.x - w.x) + (2 * k * w.x * f.x) + (2 * px.x) * ( (m * u.x) - (k * u.x) + (k * w.x) ) )/ (2 * m * w.x);
  int _py = ( (k * w.y) * (u.y - w.y) + (2 * k * w.y * f.y) + (2 * px.y) * ( (m * u.y) - (k * u.y) + (k * w.y) ) )/ (2 * m * w.y);
#ifdef DEBUG
  std::cout << "Map: ( " << _px << ", " << _py << " ) " << std::endl;  
#endif
  return T( _px, _py );
}


/** @} */ //end of group class.
