/**
 * \file feature.hpp
 *
 * \brief This file contains the prototype of feature extraction classic by level
 * and changing the code feature.
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

#ifndef FEATURE_HPP
#define FEATURE_HPP

#include <iostream> //std::cout, std::endl
#include <vector> //std::vector
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
//#include "opencv2/xfeatures2d/nonfree.hpp"
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

// Setting feature
#define _ORB_ 0
#define _KAZE_ 1
#define _SURF_ 2

/**
 * \class Feature
 *
 * \brief This class implements the Feature TAD to extract and compute features.
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
   * \fn Feature( cv::Mat img, std::vector< Level< T > > levels, int method )
   *
   * \brief Constructor default.
   *
   * \param img - Image to be foveated
   * \param levels - Fovea levels
   * \param method - Feature specification configured (see settings features)
   */
  Feature(cv::Mat img, std::vector< Level< T > > levels, int method );
    
  /**
   * \fn ~Feature()
   *
   * \brief Destructor default
   */
  ~Feature();
  
  /**
   * \fn void show( cv::Mat img, std::vector< Level< T > > levels )
   *
   * \brief Display the multiresolution structure 
   * 
   * \param img - Image to be foveated
   * \param levels - Fovea levels
   */
  void show( cv::Mat img, std::vector< Level< T > > levels );
    
private:
  //
  // Attributes
  //
  std::vector< std::vector< cv::KeyPoint > > keypoints; ///< Contains all keypoints of levels
  std::vector< cv::Mat > descriptors; ///< Contains all descriptors of levels
  std::vector< float > inliersRate; ///< Relation of inliers with position of fovea by level
  
};

#endif

/**
 * \fn Feature( cv::Mat img, std::vector< Level< T > > levels, int method )
 *
 * \brief Constructor default.
 *
 * \param img - Image to be foveated
 * \param levels - Fovea levels
 * \param method - Feature specification configured (see settings features)
 */
template <typename T, typename K>
Feature< T, K >::Feature(cv::Mat img, std::vector< Level< T > > levels, int method ){
  cv::Ptr<cv::FeatureDetector> detector;
  cv::Ptr<cv::DescriptorExtractor> descriptor;
  switch ( method ){
  case _ORB_:
    std::cout << "ORB feature actived" << std::endl;
    detector = cv::ORB::create();
    descriptor = cv::ORB::create();
    break;
  case _KAZE_:
    std::cout << "KAZE feature actived" << std::endl;
    detector = cv::KAZE::create();
    descriptor = cv::KAZE::create();
    break;
    //case _SURF_:
    //detector = cv::xfeatures2d::SURF::create(400);
    //descriptor = cv::xfeatures2d::SURF::create(400);
    //break;
  default:
    std::cout << "Feature wasn't configured" << std::endl;
    break;
  }
  std::vector< cv::KeyPoint > kp;
  cv::Mat dp, output;
  for ( int i = 0; i < levels.size(); i++ ){
    cv::Mat level = levels[i].getLevel( img );
    int64 t = cv::getTickCount();
    detector->detect ( level, kp );
    t = cv::getTickCount() - t;
    std::cout << "Feature extraction = " << t*1000/cv::getTickFrequency() << " ms, ";
    t = cv::getTickCount();
    descriptor->compute ( level, kp, dp );
    t = cv::getTickCount() - t;
    std::cout << " Feature description = " << t*1000/cv::getTickFrequency() << " ms" << std::endl;
    keypoints.push_back( kp );
    descriptors.push_back( dp );
  }
  std::vector< cv::KeyPoint >().swap( kp ); // free the memory
}
    
/**
 * \fn ~Feature()
 *
 * \brief Destructor default
 */
template <typename T, typename K>
Feature< T, K >::~Feature(){
  // free the memory
  std::vector< std::vector< cv::KeyPoint > >().swap( keypoints );
  std::vector< float >().swap( inliersRate );
  std::vector< cv::Mat >().swap( descriptors );
}
  
/**
 * \fn void show( cv::Mat img, std::vector< Level< T > > levels )
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * 
 * \param img - Image to be foveated
 * \param levels - Fovea levels
 */
template <typename T, typename K>
void 
Feature< T, K >::show( cv::Mat img, std::vector< Level< T > > levels ){
  cv::Mat output;
  /*for ( int i = 0; i < levels.size(); i++ ){
    cv::Mat level = levels[i].getLevel( img );
    cv::drawKeypoints( level, keypoints[i], output, cv::Scalar::all(-1), cv::DrawMatchesFlags::DEFAULT );
    cv::imshow( "keypoints", output);
    cv::waitKey( 0 );
  }*/
  for (int i = levels.size() - 1; i < levels.size(); i++ ){
    cv::Mat level = levels[i].getLevel( img );
    cv::drawKeypoints( level, keypoints[i], output, cv::Scalar::all(-1), cv::DrawMatchesFlags::DEFAULT );
    cv::imshow( "keypoints", output);
  }
}  
  

/** @} */ //end of group class.
