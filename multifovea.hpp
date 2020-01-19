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

using namespace std;
using namespace cv;

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
   * \fn Multifovea(Mat img, int m, T w, vector< T > fs, int mode)
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
  Multifovea(Mat img, int m, T w, vector< T > fs, int mode);
  
  /**
   * \fn Multifovea(Mat img, String ymlFile, int mode)
   *
   * \brief Constructor default of multifovea class.
   * This constructor is used to configure multiple foveas using
   * a file yaml.
   * 
   * \param img - Image to be foveated
   * \param ymlFile - File that contains all information of configuration
   * \param mode - identify the approach: reexecution, pixel-by-pixel, 
   * bitmap or block-based (see settings multifovea approaches )
   */
  Multifovea(Mat img, String ymlFile, int mode);
  
  /**
   * \fn Multifovea(int m, T w, T u, vector< T > fs, int mode)
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
  Multifovea(int m, T w, T u, vector< T > fs, int mode);
  
  /**
   * \fn ~Multifovea()
   *
   * \brief Destructor default of multifovea class.
   */
  ~Multifovea();

  /**
   * \fn void updateMultifovea(vector< T > fs)
   *
   * \brief This method update the foveas structure
   * and can be thought in update only what need to 
   * be changed, in other words, the update occors
   * only in shape of level. 
   * It is necessary to be implemented!
   *
   * \param fs - Vector with positions (x, y) of the foveas
   */
  void updateMultifovea(vector< T > fs);

  /**
   * \fn void foveatedFeatures( Mat img, int feature, int code )
   *
   * \brief This method compute and extract features
   * of foveated structures using MRMF or MMF.
   *
   * \param img - Image to be foveated
   * \param feature - This feature indicates which 
   * method feature will be choose to be extracted. 
   * ( see feature.hpp in setting feature )
   * \param code - This code indicates which method
   * to foveation will be used. If code is zero, then
   * MRMF is chosen, otherwise MMF. ( see fovea.hp )
   */
  void foveatedFeatures( Mat img, int feature, int code );
  
  /**
   * \fn Mat multifoveatedImage( Mat img )
   *
   * \brief This function builds an image with multiples focus.
   *
   * \param img - Image to be multifoveated
   *
   * \return Image multifoveated created by multiples focus
   */
  Mat multifoveatedImage( Mat img );
  
  /**
   * \fn Mat multifoveaLevelsImage( Mat img, vector< Scalar > colors )
   *
   * \brief This function builds multiples images foveated to different focus
   *
   * \param img - Image to be multifoveated
   * \param colors - Colors to paint levels
   *
   * \return Show the extraction feature in each fovea
   */
  Mat multifoveaLevelsImage( Mat img, vector< Scalar > colors );

  /**
   * \fn vector< Fovea< T >* > getFoveas()
   *
   * \brief This function return all foveas built.
   *
   * \return Vector with all foveas
   */
  vector< Fovea< T >* > getFoveas();

  /**
   * \fn Fovea< T >* getFovea( int k )
   *
   * \brief  This function return the fovea solicited
   *
   * \param k - Level of fovea
   *
   * \return The fovea of k index
   */
  Fovea< T >* getFovea( int k );
  
  /**
   * \fn void matching( vector< KeyPoint > modelKeypoints, Mat modelDescriptors )
   *
   * \brief This method realize the match between two foveas.
   *
   * \param modelKeypoints - model keypoints
   * \param modelDescriptors - model descriptors
   */
  void matching( Mat scene, Mat model, vector< KeyPoint > modelKeypoints, Mat modelDescriptors );
  
private:
  //
  // Methods
  //
  
  //
  // Attributes
  //
  vector< Fovea< T >* > foveas; ///< Pointers to foveas
  // Parameters
  int m; ///< Number levels of fovea
};

#endif

/**
 * \fn Multifovea(Mat img, int m, T w, vector< T > fs, int mode)
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
Multifovea< T >::Multifovea(Mat img, int m, T w, vector< T > fs, int mode){
  // Keeping value of quantity levels
  this->m = m;
  Fovea< T > *fovea;
  switch ( mode ){
  case REEXECUTION: {
    cout << "reexecution approach" << endl;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
    for ( int f = 0; f < fs.size(); f++ ){
      fovea = new Fovea< T >( img, m, w, fs[f] );
      foveas.push_back( fovea );
    }
    break;
  }
  case PIXELBYPIXEL: {
    cout << "pixel by pixel processing approach" << endl;
    break;
  }
  case BITMAP: {
    cout << "bitmap approach" << endl;
    break;
  }
  case BLOCKS: {
    cout << "Block-based processing approach" << endl;
    break;
  }
  default:
    cout << "There was not configured the mode" << endl;
    break;
  }
}

/**
 * \fn Multifovea(Mat img, String ymlFile, int mode)
 *
 * \brief Constructor default of multifovea class.
 * This constructor is used to configure multiple foveas using
 * a file yaml.
 * 
 * \param img - Image to be foveated
 * \param ymlFile - File that contains all information of configuration
 * \param mode - identify the approach: reexecution, pixel-by-pixel, 
 * bitmap or block-based (see settings multifovea approaches )
 */
template <typename T>
Multifovea< T >::Multifovea(Mat img, String ymlFile, int mode){
  vector< int > fx;
  FileStorage fs(ymlFile, FileStorage::READ);
  int numberOfLevels = (int) fs["numberOfLevels"];
  // Keeping value of quantity levels
  this->m = numberOfLevels - 1;
  fs["foveax"] >> fx;
  fs.release();
  Fovea< T > *fovea;
  switch ( mode ){
  case REEXECUTION: {
    cout << "reexecution approach" << endl;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
    for ( int f = 0; f < fx.size(); f++ ){
      fovea = new Fovea< T >( img, ymlFile, f );
      foveas.push_back( fovea );
    }
    break;
  }
  case PIXELBYPIXEL: {
    cout << "pixel by pixel processing approach" << endl;
    break;
  }
  case BITMAP: {
    cout << "bitmap approach" << endl;
    break;
  }
  case BLOCKS: {
    cout << "Block-based processing approach" << endl;
    break;
  }
  default:
    cout << "There was not configured the mode" << endl;
    break;
  }
}

/**
 * \fn Multifovea(int m, T w, T u, vector< T > fs, int mode)
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
Multifovea< T >::Multifovea(int m, T w, T u, vector< T > fs, int mode){
  // Keeping value of quantity levels
  this->m = m;
  Fovea< T > *fovea;
  switch ( mode ){
  case REEXECUTION: {
    cout << "reexecution approach" << endl;
#ifdef _OPENMP
#pragma omp parallel for // reference http://ppc.cs.aalto.fi/ch3/nested/
#endif
    for ( int f = 0; f < fs.size(); f++ ){
      fovea = new Fovea< T >( m, w, u, fs[f] );
      foveas.push_back( fovea );
    }
    break;
  }
  case PIXELBYPIXEL: {
    cout << "pixel by pixel processing approach" << endl;
    break;
  }
  case BITMAP: {
    cout << "bitmap approach" << endl;
    break;
  }
  case BLOCKS: {
    cout << "Block-based processing approach" << endl;
    break;
  }
  default:
    cout << "There was not configured the mode" << endl;
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
  vector< Fovea< T > >().swap( this->foveas ); // Free the memory
}

/**
 * \fn void updateMultifovea(vector< T > fs)
 *
 * \brief This method update the foveas structures
 * and can be thought in update only what need to 
 * be changed, in other words, the update occors
 * only in shape of level. 
 * It is necessary to be implemented!
 *
 * \param fs - Vector with positions (x, y) of the foveas
 */
template <typename T>
void 
Multifovea< T >::updateMultifovea( vector< T > fs){
  for ( int f = 0; f < foveas.size(); f++ )
    (this->foveas[f])->updateFovea( fs[f] );
}

/**
 * \fn void foveatedFeatures( Mat img, int feature, int code )
 *
 * \brief This method compute and extract features
 * of foveated structures using MRMF or MMF.
 *
 * \param img - Image to be foveated
 * \param feature - This feature indicates which 
 * method feature will be choose to be extracted. 
 * ( see feature.hpp in setting feature )
 * \param code - This code indicates which method
 * to foveation will be used. If code is zero, then
 * MRMF is chosen, otherwise MMF. ( see fovea.hp )
 */
template <typename T>
void
Multifovea< T >::foveatedFeatures( Mat img, int feature, int code  ){
  for ( int f = 0; f < this->foveas.size(); f++ )
    (this->foveas[f])->foveatedFeatures( img, feature, code, *(this->foveas[f]) );
}

/**
 * \fn Mat multifoveatedImage( Mat img )
 *
 * \brief This function builds an image with multiples focus.
 *
 * \param img - Image to be multifoveated
 *
 * \return Image multifoveated created by multiples focus
 */
template <typename T>
Mat 
Multifovea< T >::multifoveatedImage( Mat img ){
  vector<Scalar> colors;
  srand (time(NULL)); // Initialize random seed
  for (int i = 0; i < foveas.size(); i++){
    int r = rand() % 256; // 0 - 255
    int g = rand() % 256; // 0 - 255
    int b = rand() % 256; // 0 - 255
    colors.push_back(Scalar(b, g, r));
  }
  Mat imgMultifoveated = img.clone();
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m) // Schedule(static, m) keeps the order
#endif
  for ( int k =  0; k < this->m + 1; k++ ){ // Levels
    for ( int focus = 0; focus < this->foveas.size(); focus++ ){ // foveas
      Level< T > level = (this->foveas[focus])->getLevelFromFovea( k );
      Mat imgLevel = level.getLevel( img );
      // Mapping levels to foveated image
      vector< T > mapLevel2Image = (foveas[focus])->getMapLevel2Image( k );
      T initial = mapLevel2Image[0];
      T final = mapLevel2Image[1];
#ifdef DEBUG
      cout << "(xi, yi) = (" << initial.x << ", " << initial.y << ")" << endl;
      cout << "(xf, yf) = (" << final.x << ", " << final.y << ")" << endl;
#endif
      Rect roi = Rect( initial.x, initial.y, final.x - initial.x, final.y - initial.y );
      if ( k < m ){ // Copying levels to foveated image
	resize( imgLevel, imgLevel, Size(final.x - initial.x, final.y - initial.y), 0, 0, CV_INTER_LINEAR );
	imgLevel.copyTo( imgMultifoveated( roi ) );
      }
      else
	imgLevel.copyTo( imgMultifoveated( roi ) );
      
      // Paint rectangle in each level
      rectangle(imgMultifoveated, Point(initial.x, initial.y), Point(final.x - 1, final.y - 1), colors[focus]);
     
    }
  }
  
  return imgMultifoveated;
}


/**
 * \fn Mat multifoveaLevelsImage( Mat img, vector< Scalar > colors )
 *
 * \brief This function builds multiples images foveated to different focus
 *
 * \param img - Image to be multifoveated
 * \param colors - Colors to paint levels
 *
 * \return Show the extraction feature in each fovea
 * 
 * \note This method isn't showing +4 foveas. It's necessary modify the display!
 */
template <typename T>
Mat 
Multifovea< T >::multifoveaLevelsImage( Mat img, vector< Scalar > colors ){
  /*
    Mat a, Mat b, Mat dst // a,b loaded
    hconcat(a, b, dst) // horizontal
    vconcat(a, b, dst) // vertical
   */
  Mat imageFoveated, output, suboutput;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, this->m) // Schedule(static, m) keeps the order
#endif
  for ( int focus = 0; focus < this->foveas.size(); focus++ ){
    /*if ( focus == 0 && focus % 4 == 0 )
      suboutput = output;
    else{
      if ( focus > 0 && focus % 4 == 0 )
	vconcat( suboutput, output, suboutput );
    }*/
    
    imageFoveated = (foveas[focus])->foveatedImage( img, colors[focus] );
    if ( focus % 4 == 0 )
      output = imageFoveated;
    else
      hconcat(output, imageFoveated, output);
    
    char buffer[50];
    sprintf(buffer, "Fovea %d", focus);
    putText(output, buffer, T(10, (10*focus)+10), FONT_HERSHEY_SIMPLEX, 0.25, colors[focus], 1, 1);
    
  }
  /*if ( foveas.size() > 4 ){
    output = suboutput;
  }*/
  return output;
}

/**
 * \fn vector< Fovea< T >* > getFoveas()
 *
 * \brief This function return all foveas built.
 *
 * \return Vector with all foveas
 */
template <typename T>
vector< Fovea< T >* > 
Multifovea< T >::getFoveas(){
  return foveas;
}

/**
 * \fn Fovea< T >* getFovea( int k )
 *
 * \brief This function return the fovea solicited.
 *
 * \param k - Level of fovea
 *
 * \return The fovea of k index
 */
template <typename T>
Fovea< T >*
Multifovea< T >::getFovea( int k ){
  return foveas[k];
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
Multifovea< T >::matching( Mat scene, Mat model, vector< KeyPoint > modelKeypoints, Mat modelDescriptors ){
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for ( int focus = 0; focus < this->foveas.size(); focus++ ){
    foveas[focus]->matching( scene, model, modelKeypoints, modelDescriptors );
  }
}


/** @} */ //end of group class.
