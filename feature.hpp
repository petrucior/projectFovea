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
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"

#include "foveatedHessianDetector.hpp"

#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;


/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

// Setting feature
#define _ORB_ 0
#define _KAZE_ 1
#define _SURF_ 2
#define _AKAZE_ 3
#define _BRISK_ 4
#define _FOVEATEDFEATURES_ 5

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
   * \fn Feature( Mat img, vector< Level< T > > levels, int method )
   *
   * \brief Constructor default.
   *
   * \param img - Image to be foveated
   * \param levels - Fovea levels
   * \param method - Feature specification configured (see settings features)
   */
  Feature( Mat img, vector< Level< T > > levels, int method );

  /**
   * \fn ~Feature()
   *
   * \brief Destructor default
   */
  ~Feature();
  
  /**
   * \fn void show( Mat img, vector< Level< T > > levels )
   *
   * \brief Display the multiresolution structure 
   * 
   * \param img - Image to be foveated
   * \param levels - Fovea levels
   */
  void show( Mat img, vector< Level< T > > levels );
  
  /**
   * \fn vector< KeyPoint > getKeyPoints( int k );
   *
   * \brief This method have function to return keypoints.
   * 
   * \param k - Level of multiresolution that is solicited
   *
   * \return Return keypoints of level k requested.
   */
  vector< KeyPoint > getKeyPoints( int k );
  
  /**
   * \fn Mat getDescriptors( int k );
   *
   * \brief This method have function to return descriptors.
   * 
   * \param k - Level of multiresolution that is solicited
   *
   * \return Return descriptors of level k requested.
   */
  Mat getDescriptors( int k );
  
private:
  //
  // Attributes
  //
  vector< vector< KeyPoint > > keypoints; ///< Contains all keypoints of levels
  vector< Mat > descriptors; ///< Contains all descriptors of levels
  vector< float > inliersRate; ///< Relation of inliers with position of fovea by level  
};

#endif

/**
 * \fn Feature( Mat img, vector< Level< T > > levels, int method )
 *
 * \brief Constructor default.
 *
 * \param img - Image to be foveated
 * \param levels - Fovea levels
 * \param method - Feature specification configured (see settings features)
 */
template <typename T, typename K>
Feature< T, K >::Feature(Mat img, vector< Level< T > > levels, int method ){
  Ptr< FeatureDetector > detector;
  Ptr< DescriptorExtractor > descriptor;
  
  // ORB configuration
  int orb_nfeatures = 500;
  float orb_scaleFactor = 1.2f;
  int orb_nlevels = 1;
  int orb_edgeThreshold = 31;
  int orb_firstLevel = 0;
  int orb_WTA_K = 2;
  int orb_scoreType = ORB::HARRIS_SCORE;
  int orb_patchSize = 31;
  int orb_fastThreshold = 20;

  // KAZE Configuration
  bool kaze_extended = false;
  bool kaze_upright = false;
  float kaze_threshold = 0.001f;
  int kaze_nOctaves = 4;
  int kaze_nOctaveLayers = 4;
  int kaze_diffusivity = KAZE::DIFF_PM_G2;
  
  // AKAZE Configuration
  int akaze_descriptor_type = AKAZE::DESCRIPTOR_MLDB;
  int akaze_descriptor_size = 0;
  int akaze_descriptor_channels = 3;
  float akaze_threshold = 0.001f;
  int akaze_nOctaves = 4;
  int akaze_nOctaveLayers = 4;
  int akaze_diffusivity = KAZE::DIFF_PM_G2; 

  // BRISK Configuration
  int thresh = 30;
  int octaves = 3;
  float patternScale = 1.0f;

  // Detector
  switch ( method ){
  case _ORB_ :
#ifdef DEBUG
    cout << "ORB feature actived" << endl;
#endif
    detector = ORB::create(orb_nfeatures, orb_scaleFactor, orb_nlevels, orb_edgeThreshold, orb_firstLevel, orb_WTA_K, orb_scoreType, orb_patchSize, orb_fastThreshold);
    descriptor = ORB::create(orb_nfeatures, orb_scaleFactor, orb_nlevels, orb_edgeThreshold, orb_firstLevel, orb_WTA_K, orb_scoreType, orb_patchSize, orb_fastThreshold);
    break;
  case _KAZE_:
#ifdef DEBUG
    cout << "KAZE feature actived" << endl;
#endif
    detector = KAZE::create(kaze_extended, kaze_upright, kaze_threshold, kaze_nOctaves, kaze_nOctaveLayers, kaze_diffusivity);
    descriptor = KAZE::create(kaze_extended, kaze_upright, kaze_threshold, kaze_nOctaves, kaze_nOctaveLayers, kaze_diffusivity);
    break;
  case _SURF_:
#ifdef DEBUG
    cout << "SURF feature actived" << endl;
#endif
    detector = SURF::create(400);
    descriptor = SURF::create(400);
    break;
  case _AKAZE_:
#ifdef DEBUG
    cout << "AKAZE feature actived" << endl;
#endif
    detector = AKAZE::create(akaze_descriptor_type, akaze_descriptor_size, akaze_descriptor_channels, akaze_threshold, akaze_nOctaves, akaze_nOctaveLayers, akaze_diffusivity);
    descriptor = AKAZE::create(akaze_descriptor_type, akaze_descriptor_size, akaze_descriptor_channels, akaze_threshold, akaze_nOctaves, akaze_nOctaveLayers, akaze_diffusivity);
    break;
  case _BRISK_:
#ifdef DEBUG
    cout << "BRISK feature actived" << endl;
#endif
    detector = BRISK::create( thresh, octaves, patternScale );
    descriptor = BRISK::create( thresh, octaves, patternScale );
    break;
  case _FOVEATEDFEATURES_:
#ifdef DEBUG
    cout << "FoveatedFeatures actived" << endl;
#endif
    detector = SURF::create(400);
    descriptor = SURF::create(400);
    break;
  default:
    cout << "Feature wasn't configured" << endl;
    break;
  }
  vector< KeyPoint > kp;
  Mat dp, output;
  for ( int i = 0; i < levels.size(); i++ ){
    Mat level = levels[i].getLevel( img );
    int64 t = getTickCount();
    if ( method == _FOVEATEDFEATURES_ )
      foveatedHessianDetector(level, Mat(), kp, levels[i]);
    else
      detector->detect ( level, kp );
    t = getTickCount() - t;
#ifdef DEBUG
    cout << "Feature extraction = " << t*1000/getTickFrequency() << " ms, ";
#endif
    t = getTickCount();
    descriptor->compute ( level, kp, dp );
    t = getTickCount() - t;
#ifdef DEBUG
    cout << " Feature description = " << t*1000/getTickFrequency() << " ms" << endl;
#endif
    keypoints.push_back( kp );
    descriptors.push_back( dp );

    vector< KeyPoint >().swap( kp ); // free the memory
  }
  //vector< KeyPoint >().swap( kp ); // free the memory
}

/**
 * \fn ~Feature()
 *
 * \brief Destructor default
 */
template <typename T, typename K>
Feature< T, K >::~Feature(){
  // free the memory
  vector< vector< KeyPoint > >().swap( keypoints );
  vector< float >().swap( inliersRate );
  vector< Mat >().swap( descriptors );
}
  
/**
 * \fn void show( Mat img, vector< Level< T > > levels )
 *
 * \brief Pure virtual method caracterize this class like abstract.
 * 
 * \param img - Image to be foveated
 * \param levels - Fovea levels
 */
template <typename T, typename K>
void 
Feature< T, K >::show( Mat img, vector< Level< T > > levels ){
  Mat output;
  /*for ( int i = 0; i < levels.size(); i++ ){
    Mat level = levels[i].getLevel( img );
    drawKeypoints( level, keypoints[i], output, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    imshow( "keypoints", output);
    waitKey( 0 );
  }*/
  for (int i = levels.size() - 1; i < levels.size(); i++ ){
    Mat level = levels[i].getLevel( img );
    drawKeypoints( level, keypoints[i], output, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    imshow( "keypoints", output);
  }
}  

/**
 * \fn vector< KeyPoint > getKeyPoints( int k );
 *
 * \brief This method have function to return keypoints.
 * 
 * \param k - Level of multiresolution that is solicited
 *
 * \return Return keypoints of level k requested.
 */
template <typename T, typename K>
vector< KeyPoint > 
Feature< T, K >::getKeyPoints( int k ){
  return this->keypoints[k];
}

/**
 * \fn Mat getDescriptors( int k );
 *
 * \brief This method have function to return descriptors.
 * 
 * \param k - Level of multiresolution that is solicited
 *
 * \return Return descriptors of level k requested.
 */
template <typename T, typename K>
Mat
Feature< T, K >::getDescriptors( int k ){
  return this->descriptors[k];
}

/** @} */ //end of group class.
