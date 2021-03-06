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
#include <stdio.h>
#include <vector> //std::vector
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"

#include "level.hpp" //std::vector< Level >
#include "feature.hpp"

#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

// Settings fovea
#define MRMF 0 ///< Identify the use of MRMF method
#define MMF 1 ///< Identify the use of MMF method

using namespace std;
using namespace cv;

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
   * \fn Fovea( int m, T w, T u, T f, int shapeMode )
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
   * \param shapeMode - Feature specification configured (see settings 
   * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
   */
  Fovea( int m, T w, T u, T f, int shapeMode );

  /**
   * \fn Fovea( T u, String ymlFile, int index, int shapeMode )
   *
   * \brief Constructor default of fovea class.
   * This constructor is used to configure multiple foveas using
   * a file yaml.
   *
   * \param u - Size of image
   * \param ymlFile - File that contains all information of configuration
   * \param index - Vector index with fovea position information
   * \param shapeMode - Feature specification configured (see settings 
   * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
   */
  Fovea( T u, String ymlFile, int index, int shapeMode );
  
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
   * \fn const vector< T > getParameters();
   *
   * \brief Function that is responsible for informing the user of the fovea parameters
   *
   * \return The parameters of the fovea
   *
   * \note The parameter "m" is replicated to coordinates of type T
   */
  const vector< T > getParameters();

  /**
   * \fn const vector< vector< int > > getVectorsFeature();
   *
   * \brief Function that is responsible for informing bvector, etavector and levelvector
   *
   * \return The parameters bvector, etavector and levelvector of the feature
   */
  const vector< vector< int > > getVectorsFeature();
  
  /**
   * \fn const vector< double > getParametersFeature();
   *
   * \brief Function that is responsible for informing nOctaveLayers and hessianThreshold
   *
   * \return The parameters nOctaveLayers and hessianThreshold of the feature
   */
  const vector< double > getParametersFeature();
  
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
   * \fn Mat foveatedImage( Mat img, Scalar color )
   *
   * \brief This function builds the focused image.
   *
   * \param img - Image to be foveated
   * \param color - Color to paint levels
   *
   * \return Image foveated created by levels
   */
  Mat foveatedImage( Mat img, Scalar color );

  /**
   * \fn Mat foveatedImage( Mat img, Scalar color, int code )
   *
   * \brief This function builds the focused image.
   *
   * \param img - Image to be foveated
   * \param color - Color to paint levels
   * \param code - This code indicates which method
   * to foveation will be used. If code is zero, then
   * MRMF is chosen, otherwise MMF. ( see fovea.hp )
   *
   * \return Image foveated created by levels
   */
  Mat foveatedImage( Mat img, Scalar color, int code );
  
  /**
   * \fn void foveatedFeatures( Mat img, int feature, int code, Fovea< T > fovea )
   *
   * \brief This method compute and extract features 
   * of foveated structure using MRMF or MMF.
   *
   * \param img - Image to be foveated
   * \param feature - This feature indicates which 
   * method feature will be choose to be extracted. 
   * ( see feature.hpp in setting feature )
   * \param code - This code indicates which method
   * to foveation will be used. If code is zero, then
   * MRMF is chosen, otherwise MMF. ( see fovea.hp )
   * \param fovea - Fovea pointer
   */
  void foveatedFeatures( Mat img, int feature, int code, Fovea< T > fovea );
  
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
   * \fn vector< T > getMapLevel2Image( int k )
   *
   * \brief This method calculates the position of pixel on the 
   * level to image by level
   *
   * \param k - level of fovea
   *
   * \return Vector containing the boundingBox to map of level
   */
  vector< T > getMapLevel2Image( int k );
  
  /**
   * \fn Feature< T, Fovea< T > >* getFeatures()
   *
   * \brief This method return pointer to features
   *
   * \return Pointer from Features
   */
  Feature< T, Fovea< T > >* getFeatures();
  
  /**
   * \fn void matching( vector< KeyPoint > modelKeypoints, Mat modelDescriptors )
   *
   * \brief This method realize the match between two foveas.
   *
   * \param modelKeypoints - model keypoints
   * \param modelDescriptors - model descriptors
   */
  void matching( Mat scene, Mat model, vector< KeyPoint > modelKeypoints, Mat modelDescriptors );
  
  /**
   * \fn double getInliersRatio( int k )
   *
   * \brief This method returns inlier ratio
   *
   * \param k - Level from fovea
   *
   * \return This method return inlier ratio calculated after features extraction and match in level k
   */
  double getInliersRatio( int k );
  
  /**
   * \fn int getNumberMatches( int k )
   *
   * \brief This method returns number of matches
   *
   * \param k - Level from fovea
   *
   * \return This method return the number of matches calculated after features extraction and match in level k
   */
  int getNumberMatches( int k );
  
  
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
  vector< Level< T > > levels; ///< List of levels
  Feature< T, Fovea< T > >* features = NULL; ///< Features
  vector< double > inliersRatio; ///< Inliers Ratio
  vector< int > numberMatches; ///< Quantity of matches
  int m; ///< Number levels of fovea
  T w; ///< Size of levels
  T u; ///< Size of image
  T f; ///< Position (x, y) to build the fovea
  vector< int > bvector; // bvector: [b1, b2, ..., bn]: a vector where bi is 0 if the feature extraction step number i should be discarded or 1, otherwise
  vector< int > etavector; // etavector: [e1, e2, ..., en]: a vector where ei is the octave (> 0) for which the feature extraction step number i should be performed
  vector< int > levelvector; // levelvector: [l1, l2, ..., ln]: a vector where li is the foveated model level (>= 0 and < numberOfLevels) for which the feature extraction step number should be performed
  int nOctaveLayers;
  int hessianThreshold;
  int shape = 0; ///< Initialization shape
};

#endif

/**
 * \fn Fovea( int m, T w, T u, T f, int shapeMode )
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
 * \param shapeMode - Feature specification configured (see settings 
 * in level.hpp), where 0 (blocks), 1 (rectangle) or 2 (polygons)
 */
template <typename T>
Fovea< T >::Fovea( int m, T w, T u, T f, int shapeMode ){
  // Cleaning parameters
  bvector.clear();
  etavector.clear();
  levelvector.clear();
  nOctaveLayers = 3;
  hessianThreshold = 100;
  shape = shapeMode;
  
  this->checkParameters( m, w, u, f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m + 1; k++ ){
    Level< T > l( k, m, w, u, f, shapeMode );
    levels.push_back( l );
  }
}

/**
 * \fn Fovea( T u, String ymlFile, int index, int shapeMode )
 *
 * \brief Constructor default of fovea class.
 * This constructor is used to configure multiple foveas using
 * a file yaml.
 *
 * \param u - Size of image
 * \param ymlFile - File that contains all information of configuration
 * \param index - Vector index with fovea position information
 * \param shapeMode - Feature specification configured (see settings 
 * in level.hpp), where 0 (blocks), 1 (rectangle) or 2 (polygons)
 */
template <typename T>
Fovea< T >::Fovea( T u, String ymlFile, int index, int shapeMode ){
  // Cleaning parameters
  bvector.clear();
  etavector.clear();
  levelvector.clear();
  nOctaveLayers = 3;
  hessianThreshold = 100;
  shape = shapeMode;
    
  FileStorage fs(ymlFile, FileStorage::READ);
  //int ux = (int) fs["imageWidth"]; // img.cols
  //int uy = (int) fs["imageHeight"]; // img.rows
  //T u = T( ux, uy );
  int wx = (int) fs["smallestLevelWidth"];
  int wy = (int) fs["smallestLevelHeight"];
  T w = T( wx, wy );
  fs["bvector"] >> bvector;
  fs["etavector"] >> etavector;
  fs["levelvector"] >> levelvector;
  int numberOfLevels = (int) fs["numberOfLevels"];
  m = numberOfLevels - 1;
  fs["hessianThreshold"] >> hessianThreshold;
  fs["nOctaveLayers"] >> nOctaveLayers;
  vector< int > fx, fy;
  fs["foveax"] >> fx;
  fs["foveay"] >> fy;
  // handling possible problems with foveas vector index
  if ( fx.size() < index + 1 )
    return;
  T f = T( fx[index], fy[index] );
  fs.release();
  this->checkParameters( m, w, u, f );
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
  for ( int k = 0; k < m + 1; k++ ){
    Level< T > l( k, m, w, u, f, shapeMode );
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
  vector< Level< T > >().swap( this->levels ); // Free the memory
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
 * \fn const vector< T > getParameters();
 *
 * \brief Function that is responsible for informing the user of the fovea parameters
 *
 * \return The parameters of the fovea
 *
 * \note The parameter "m" is replicated to coordinates of type T
 */
template <typename T>
const vector< T > 
Fovea< T >::getParameters(){
  vector< T > parameters;
  parameters.push_back( T( this->m, this->m ) );
  parameters.push_back( this->w );
  parameters.push_back( this->u );
  parameters.push_back( this->f );
  return parameters;
}

/**
 * \fn const vector< vector< int > > getVectorsFeature();
 *
 * \brief Function that is responsible for informing bvector, etavector and levelvector
 *
 * \return The parameters bvector, etavector and levelvector of the feature
 */
template <typename T>
const vector< vector< int > >
Fovea< T >::getVectorsFeature(){
  vector< vector< int > > parameters;
  parameters.push_back( bvector );
  parameters.push_back( etavector );
  parameters.push_back( levelvector );
  return parameters;
}

/**
 * \fn const vector< double > getParametersFeature();
 *
 * \brief Function that is responsible for informing nOctaveLayers and hessianThreshold
 *
 * \return The parameters nOctaveLayers and hessianThreshold of the feature
 */
template <typename T>
const vector< double >
Fovea< T >::getParametersFeature(){
  vector< double > parameters;
  parameters.push_back( nOctaveLayers );
  parameters.push_back( hessianThreshold );
  return parameters;
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
  for ( int k = 0; k < m + 1; k++ ){
    levels[k].updateLevel( m, w, u, this->f, shape );
  }
}

/**
 * \fn Mat foveatedImage( Mat img, Scalar color )
 *
 * \brief This function builds the focused image.
 *
 * \param img - Image to be foveated
 * \param color - Color to paint levels
 *
 * \return Image foveated created by levels
 */
template <typename T>
Mat 
Fovea< T >::foveatedImage( Mat img, Scalar color ){
  Mat imgFoveated = img.clone();
  vector< KeyPoint > kp;
  vector< vector< KeyPoint > > keypoints;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m) // Schedule(static, m) keeps the order
#endif
  for ( int k = 0; k < levels.size(); k++ ){ // Levels
    Mat imgLevel = levels[k].getLevel( img );
    if ( features != NULL ){
      kp = features->getKeyPoints( k );
      for ( int i = 0; i < kp.size(); i++ ){
	Point2f kpPos = mapLevel2Image( k, this->m, this->w, this->u, this->f, Point2f( kp[i].pt.x, kp[i].pt.y ) );
	kp[i].pt.x = kpPos.x;
	kp[i].pt.y = kpPos.y;
      }
      keypoints.push_back( kp );
    }
    // Mapping levels to foveated image
    T initial = mapLevel2Image( k, this->m, this->w, this->u, this->f, T( 0, 0 ) ); 
    T final = mapLevel2Image( k, this->m, this->w, this->u, this->f, T( this->w.x, this->w.y ) );
#ifdef DEBUG
    cout << "(xi, yi) = (" << initial.x << ", " << initial.y << ")" << endl;
    cout << "(xf, yf) = (" << final.x << ", " << final.y << ")" << endl;
#endif
    Rect roi = Rect( initial.x, initial.y, final.x - initial.x, final.y - initial.y );
    if ( k < m ){ // Copying levels to foveated image
      resize( imgLevel, imgLevel, Size(final.x - initial.x, final.y - initial.y), 0, 0, CV_INTER_LINEAR );
      //imgLevel.copyTo( imgFoveated( roi ) );
    }
    //else
      imgLevel.copyTo( imgFoveated( roi ) );
    
    // Paint rectangle in each level
    rectangle(imgFoveated, Point(initial.x, initial.y), Point(final.x - 1, final.y - 1), color);
    
  }
  
  if ( keypoints.size() != 0 ){
    for ( int i = 0; i < keypoints.size(); i++ )
      drawKeypoints( imgFoveated, keypoints[i], imgFoveated, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
  }
  
  features = NULL;
  return imgFoveated;
}

/**
 * \fn Mat foveatedImage( Mat img, Scalar color, int code )
 *
 * \brief This function builds the focused image.
 *
 * \param img - Image to be foveated
 * \param color - Color to paint levels
 * \param code - This code indicates which method
 * to foveation will be used. If code is zero, then
 * MRMF is chosen, otherwise MMF. ( see fovea.hp )
 *
 * \return Image foveated created by levels
 */
template <typename T>
Mat
Fovea< T >::foveatedImage( Mat img, Scalar color, int code ){
  Mat imgFoveated = img.clone();
  if ( code == MRMF ){
    //cout << "MRMF actived" << endl;
    imgFoveated = foveatedImage( img, color );
  }
  if ( code == MMF ){
    //cout << "MMF actived" << endl;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m) // Schedule(static, m) keeps the order
#endif
    for ( int k = 0; k < levels.size(); k++ ){ // Levels
      vector< T > params = levels[k].boundingBox();
      int dx = params[0].x;
      int dy = params[0].y;
      int sx = params[1].x;
      int sy = params[1].y;
      rectangle( imgFoveated, Point(dx, dy), Point(dx+sx, dy+sy), Scalar(255, 255, 255) );
    }
    vector< KeyPoint > kp = features->getKeyPoints( 0 );
    if ( kp.size() != 0 )
      drawKeypoints( imgFoveated, kp, imgFoveated, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
  }
  return imgFoveated;
}

/**
 * \fn bool foveatedFeatures( Mat img, int feature, int code, Fovea< T >* fovea )
 *
 * \brief This method compute and extract features 
 * of foveated structure using MRMF or MMF.
 *
 * \param img - Image to be foveated
 * \param feature - This feature indicates which 
 * method feature will be choose to be extracted. 
 * ( see feature.hpp in setting feature )
 * \param code - This code indicates which method
 * to foveation will be used. If code is zero, then
 * MRMF is chosen, otherwise MMF. ( see fovea.hp )
 * \param fovea - Fovea pointer
 */
template <typename T>
void
Fovea< T >::foveatedFeatures( Mat img, int feature, int code, Fovea< T > fovea ){
  if ( code == MRMF ){
    //cout << "MRMF actived" << endl;
    features = new Feature< T, Fovea< T > >( img, levels, feature );
#ifdef DEBUG
    features->show( img, levels );
#endif
  }
  if ( code == MMF ){
    //cout << "MMF actived" << endl;
    features = new Feature< T, Fovea< T > >( img, fovea, feature );
  }
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
 * \fn vector< T > getMapLevel2Image( int k )
 *
 * \brief This method calculates the position of pixel on the 
 * level to image by level
 *
 * \param k - level of fovea
 *
 * \return Vector containing the boundingBox to map of level
 */
template <typename T>
vector< T > 
Fovea< T >::getMapLevel2Image( int k ){
  vector< T > mapImage;
  mapImage.push_back( this->mapLevel2Image( k, this->m, this->w, this->u, this->f, T( 0, 0 ) ) );
  mapImage.push_back( this->mapLevel2Image( k, this->m, this->w, this->u, this->f, T( this->w.x, this->w.y ) ) );
  return mapImage;
}

/**
 * \fn Feature< T, Fovea< T > >* getFeatures()
 *
 * \brief This method return pointer to features
 *
 * \return Pointer from Features
 */
template <typename T>
Feature< T, Fovea< T > >* 
Fovea< T >::getFeatures(){
  return features;
}

/**
 * \fn void matching( vector< KeyPoint > modelKeypoints, Mat modelDescriptors )
 *
 * \brief This method realize the match between two foveas.
 *
 * \param modelKeypoints - model keypoints
 * \param modelDescriptors - model descriptors
 */
template <typename T>
void
Fovea< T >::matching( Mat scene, Mat model, vector< KeyPoint > modelKeypoints, Mat modelDescriptors ){
  // -----------------------
  // First possibility
  // -----------------------
  /*Ptr<BFMatcher> matcher = BFMatcher::create ( NORM_HAMMING, true );
  vector< DMatch > matches;
  vector< Point2f > modelPoints, imgPoints;
  Mat mask, level;
  int inliers, outliers;
  inliersRatio.clear();
  numberMatches.clear();
  //int64 t = getTickCount();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->levels.size()) // Schedule(static, m+1) keeps the order
#endif
  for ( int k = 0; k < levels.size(); k++ ){ // Levels
    inliers = 0; outliers = 0;
#ifdef DEBUG
    cout << "level " << k << endl;
#endif
    if ( !(this->features)->getDescriptors( k ).empty() ){
      level = levels[k].getLevel( scene );
      matches.clear();
      matcher->match( modelDescriptors, (this->features)->getDescriptors( k ), matches );
      
      //Mat img_matches;
      //drawMatches(model, modelKeypoints, level, (this->features)->getKeyPoints( k ), matches, img_matches);
      //imshow( "img_matches", img_matches );
      //waitKey( 0 );

      if ( matches.size() != 0 ){
	// -------------------
	// First tentative
	// -------------------
	//modelPoints.clear(); imgPoints.clear();
	//vector< KeyPoint > imgKeypoints =  (this->features)->getKeyPoints( k );
	//for ( unsigned int i = 0; i < matches.size(); i++ ){
	//  DMatch m = matches[i];
	//  modelPoints.push_back(modelKeypoints[m.queryIdx].pt);
	//  imgPoints.push_back(imgKeypoints[m.trainIdx].pt);
	//}
	
	//
	// --------------------------------------------------
	// This function is very expensive to be processing
	// --------------------------------------------------
	//
	//Mat H = findHomography( modelPoints, imgPoints, RANSAC, 2.5f, mask );
	////Mat H = findHomography( modelPoints, imgPoints, mask, RANSAC, 3 );
	//for ( int x = 0; x < mask.rows; x++ ){
	//  for ( int y = 0; y < mask.cols; y++ ){
	//    //Counting inliers and outliers
	//    (int)mask.at<uchar>(x, y) == 1 ? inliers++ : outliers++;
	//  }
	//}
	
	
	// -------------------
	// Second tentative
	// -------------------
	// Detectin min and max distances between matches
	//double min_dist=matches[0].distance, max_dist=matches[0].distance;
	//for ( int i = 1; i < modelDescriptors.rows; i++ ){
	//  double dist = matches[i].distance;
	//  if ( dist < min_dist ) min_dist = dist;
	//  if ( dist > max_dist ) max_dist = dist;
	//}
	
	// When the distance between descriptors is greater than 2x the minimum distance, 
	// the match is considered incorrect, but sometimes the minimum distance is too 
	// small and an empirical value of 30 is set to the lower limit.
	//vector< DMatch > good_matches;
	//for ( int i = 0; i < modelDescriptors.rows; i++ ){
	//  if ( matches[i].distance <= max( 2*min_dist, 30.0 ) ){
	//    good_matches.push_back ( matches[i] );
	//    inliers++;
	//  }
	//  else
	//    outliers++;
	//}

	// -------------------
	// Third tentative
	// -------------------
	double min_dist=matches[0].distance, max_dist=matches[0].distance;
	for ( int i = 1; i < modelDescriptors.rows; i++ ){
	  double dist = matches[i].distance;
	  if ( dist < min_dist ) min_dist = dist;
	  if ( dist > max_dist ) max_dist = dist;
	}
	
	modelPoints.clear(); imgPoints.clear();
	vector< KeyPoint > imgKeypoints =  (this->features)->getKeyPoints( k );
	for ( unsigned int i = 0; i < matches.size(); i++ ){
	  if ( matches[i].distance <= max( 2*min_dist, 30.0 ) ){
	    DMatch m = matches[i];
	    modelPoints.push_back(modelKeypoints[m.queryIdx].pt);
	    imgPoints.push_back(imgKeypoints[m.trainIdx].pt);
	  }
	}
	
	if ( modelPoints.size() != 0 ){
	  Mat H = findHomography( modelPoints, imgPoints, RANSAC, 2.5f, mask );
	  for ( int x = 0; x < mask.rows; x++ ){
	    for ( int y = 0; y < mask.cols; y++ ){
	      //Counting inliers and outliers
	      (int)mask.at<uchar>(x, y) == 1 ? inliers++ : outliers++;
	    }
	  }
	}
	
#ifdef DEBUG
	cout << "quantidade de inliers: " << inliers << endl;
	cout << "quantidade de outliers: " << outliers << endl;
#endif
      }
    }
    else{ // Descriptors empty
#ifdef DEBUG
      cout << "quantidade de inliers: " << inliers << endl;
      cout << "quantidade de outliers: " << outliers << endl;
#endif
    }
    
    double result = 0.0;
    if ( matches.size() != 0 )
      result = (inliers * 1.0)/matches.size();
    inliersRatio.push_back( result );
    numberMatches.push_back( matches.size() );
    
  }

  //t = getTickCount() - t;
  //cout << "time = " << t*1000/getTickFrequency() << " ms ";
  */
  
  
  // -----------------------
  // Second possibility
  // -----------------------
  /*
  const double ransac_thresh = 2.5f; // RANSAC inlier threshold
  const double nn_match_ratio = 0.8f; // Nearest-neighbour matching ratio
  Ptr< DescriptorMatcher> matcher  = DescriptorMatcher::create( "BruteForce-Hamming" );
  vector< vector< DMatch > > matches;
  vector< Point2f > modelPoints, imgPoints;
  Mat mask, level;
  int inliers, outliers;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->levels.size()) // Schedule(static, m+1) keeps the order
#endif
  //for ( int k = 0; k < 1; k++ ){ // Levels
  for ( int k = 0; k < levels.size(); k++ ){ // Levels
    inliers = 0; outliers = 0;
    //#ifdef DEBUG
    cout << "level " << k << endl;
    //#endif
    if ( !(this->features)->getDescriptors( k ).empty() ){
      level = levels[k].getLevel( scene );
      matches.clear();
      matcher->knnMatch( modelDescriptors, (this->features)->getDescriptors( k ), matches, 20 );
      
      if ( matches.size() != 0 ){
	modelPoints.clear(); imgPoints.clear();
	vector< KeyPoint > imgKeypoints =  (this->features)->getKeyPoints( k );
	for ( unsigned int i = 0; i < matches.size(); i++ ){
	  if( matches[i][0].distance < nn_match_ratio * matches[i][1].distance ) {
	    DMatch m = matches[i][0];
	    modelPoints.push_back(modelKeypoints[m.queryIdx].pt);
	    imgPoints.push_back(imgKeypoints[m.trainIdx].pt);
	  }
	}
	if ( modelPoints.size() >= 4 ){
	  Mat H = findHomography( modelPoints, imgPoints, RANSAC, ransac_thresh, mask);
	  //cout << " , mask " << mask.rows << ", " << mask.cols << endl;
	  for ( int x = 0; x < mask.rows; x++ ){
	    for ( int y = 0; y < mask.cols; y++ ){
	      //cout << (int)mask.at<uchar>(x,y) << endl;
	      //Counting inliers and outliers
	      (int)mask.at<uchar>(x, y) == 1 ? inliers++ : outliers++;
	    }
	  }
	}
	//#ifdef DEBUG
	cout << "quantidade de inliers: " << inliers << endl;
	cout << "quantidade de outliers: " << outliers << endl;
	//#endif
      }
    }
    else{
      //#ifdef DEBUG
      cout << "quantidade de inliers: " << inliers << endl;
      cout << "quantidade de outliers: " << outliers << endl;
      //#endif
    }
    double matchesTotal = inliers + outliers;
    double result = 0.0;
    if ( matchesTotal > 0.0 )
      result = ((inliers * 1.0)/matchesTotal);
    inliersRatio.push_back( result );
    numberMatches.push_back( matches.size() );
  }
  */


  // -----------------------
  // Third possibility
  // -----------------------
  //BFMatcher matcher( NORM_HAMMING, true); // https://docs.opencv.org/3.4/d3/da1/classcv_1_1BFMatcher.html
  //BFMatcher matcher( NORM_L2, true );
  //Ptr< BFMatcher > matcher = BFMatcher::create( NORM_L2, true );
  Ptr< DescriptorMatcher > matcher = DescriptorMatcher::create( DescriptorMatcher::FLANNBASED );
  vector< DMatch > matches;
  vector< Point2f > modelPoints, imgPoints;
  Mat mask, level;
  int inliers, outliers;
  int matchesMax = 0;
  inliersRatio.clear();
  numberMatches.clear();
  vector< DMatch > good_matches;
  string type = "bf"; //"bf/knn";
  //int64 t = getTickCount();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->levels.size()) // Schedule(static, m+1) keeps the order
#endif
  for ( int k = 0; k < levels.size(); k++ ){ // Levels
    inliers = 0; outliers = 0;
#ifdef DEBUG
    cout << "level " << k << endl;
#endif
    good_matches.clear();
    if ( !(this->features)->getDescriptors( k ).empty() ){
      level = levels[k].getLevel( scene );
      matches.clear();
      
      vector< KeyPoint > sceneKeypoints =  (this->features)->getKeyPoints( k );
      if ( sceneKeypoints.size() > 2 ){
	if (type == "bf") {
	  matcher->match(modelDescriptors, (this->features)->getDescriptors( k ), matches, Mat());
	  //-- Filter matches using the Lowe's ratio test
	  const float ratio_thresh = 0.25f;
      	  for (int i = 0; i < static_cast<int>(matches.size()); ++i){
	    if ( matches[i].distance < ratio_thresh )
              good_matches.push_back(matches[i]);
	  }
	}
	if (type == "knn") {
	  vector< vector<DMatch> > vmatches;
	  matcher->knnMatch( modelDescriptors, (this->features)->getDescriptors( k ), vmatches, 2 );
	  //-- Filter matches using the Lowe's ratio test
	  const float ratio_thresh = 0.75f;
      	  for (int i = 0; i < static_cast<int>(vmatches.size()); ++i) {
	    if (!vmatches[i].size()) {
	      continue;
	    }	    
	    if (vmatches[i][0].distance < ratio_thresh * vmatches[i][1].distance)
	      good_matches.push_back(vmatches[i][0]);
	  }
	}
	
	//Mat img_matches;
	//drawMatches(model, modelKeypoints, level, (this->features)->getKeyPoints( k ), matches, img_matches);
	//imshow( "img_matches", img_matches );
	//waitKey( 0 );

	//
	// Used when you need more precision
	//
	/*sort(good_matches.begin(), good_matches.end());
	
	modelPoints.clear(); imgPoints.clear();
	for ( unsigned int i = 0; i < good_matches.size(); i++ ){
	  DMatch m = good_matches[i];
	  modelPoints.push_back(modelKeypoints[m.queryIdx].pt);
	  imgPoints.push_back(sceneKeypoints[m.trainIdx].pt);
	}
	
	if ( modelPoints.size() > 5 ){
	  Mat H = findHomography( modelPoints, imgPoints, RANSAC, 2.5f, mask);
	  for ( int x = 0; x < mask.rows; x++ ){
	    for ( int y = 0; y < mask.cols; y++ ){
	      //Counting inliers and outliers
	      (int)mask.at<uchar>(x, y) == 1 ? inliers++ : outliers++;
	    }
	  }
	}*/
	
      }
    }
#ifdef DEBUG
    cout << "quantidade de features: " << good_matches.size() << endl;
    cout << "quantidade de inliers: " << inliers << endl;
    cout << "quantidade de outliers: " << outliers << endl;
#endif
    
    /*double result = 0.0;
    if ( good_matches.size() != 0 )
      result = (double)(inliers * 1.0)/good_matches.size();
    inliersRatio.push_back( result );
    numberMatches.push_back( good_matches.size() );*/

    numberMatches.push_back( good_matches.size() );
    matchesMax += static_cast<int>(good_matches.size());
  }
  
  for( int i = 0; i < numberMatches.size(); i++ )
    inliersRatio.push_back( (double)(numberMatches[i] * 1.0 )/matchesMax );

  /*
   *  matchesMax ------ 100
   *  numberMatches[i] ----- x
   *  x = numberMatches[i] * 100 / matchesMax
   */			   
			   
  //t = getTickCount() - t;
  //cout << "time = " << t*1000/getTickFrequency() << " ms ";
  
}

/**
 * \fn double getInliersRatio( int k )
 *
 * \brief This method returns inlier ratio
 *
 * \param k - Level from fovea
 *
 * \return This method return inlier ratio calculated after features extraction and match in level k
 */
template <typename T>
double 
Fovea< T >::getInliersRatio( int k ){
  if ( inliersRatio.size() == 0 ) 
    return 0.0;
  return inliersRatio[k];
}

/**
 * \fn int getNumberMatches( int k )
 *
 * \brief This method returns number of matches
 *
 * \param k - Level from fovea
 *
 * \return This method return the number of matches calculated after features extraction and match in level k
 */
template <typename T>
int 
Fovea< T >::getNumberMatches( int k ){
  if ( numberMatches.size() == 0 )
    return 0;
  return numberMatches[k];
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
  cout << "Map: ( " << _px << ", " << _py << " ) " << endl;  
#endif
  return T( _px, _py );
}


/** @} */ //end of group class.
