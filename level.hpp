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
#include "opencv2/imgproc/imgproc.hpp"
#include "shape.hpp" //Shape
#include "rectangle.hpp" // Shape< T >::Rectangle
#include "polygons.hpp" // Shape< T >::Polygons
#include "blocks.hpp" // Shape< T >::Blocks



using namespace std;
using namespace cv;

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

// Setting Shape
#define _BLOCKS_ 0
#define _RECTANGLE_ 1
#define _POLYGONS_ 2

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
   * \fn Level( int k, int m, T w, T u, T f, int shapeMode )
   *
   * \brief Default constructor for Level class which will create 
   * a level of anyone shape.
   *
   * \param k - Level of fovea
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) to build the fovea
   * \param shapeMode - Feature specification configured (see settings 
   * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
   */
  Level( int k, int m, T w, T u, T f, int shapeMode );
  
  /**
   * \fn ~Level()
   *
   * \brief Default destructor class
   */
  ~Level();

  /**
   * \fn void updateLevel( int m, T w, T u, T f, int shapeMode )
   *
   * \brief Default constructor for Level class which will create 
   * a level of anyone shape.
   *
   * \param m - Number levels of fovea
   * \param w - Size of levels
   * \param u - Size of image
   * \param f - Position (x, y) to build the fovea
   * \param shapeMode - Feature specification configured (see settings 
   * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
   */
  void updateLevel( int m, T w, T u, T f, int shapeMode );

  /**
   * \fn void updateLevel( vector< Level< T > > levels  )
   *
   * \brief Default constructor for Level class which will create 
   * a level of Blocks shape.
   *
   * \param levels - List of all previously processed levels
   */
  void updateLevel( vector< Level< T > > levels );
  
  /**
   * \fn Mat getLevel( Mat img )
   *
   * \brief Method responsable to create an image with 
   * with level dimension
   *
   * \param img - Image will be foveated
   *
   * \return Image that represent the level
   */
  Mat getLevel( Mat img );

  /**
   * \fn vector< Mat > getLevelParts( Mat img )
   *
   * \brief Method responsable to create an image with 
   * with parts of level
   *
   * \param img - Image will be foveated
   *
   * \return Image that represent the level
   */
  vector< Mat > getLevelParts( Mat img );
  
  /**
   * \fn vector< T > boundingBox( int k, int m, T w, T u, T f )
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
  vector< T > boundingBox( int k, int m, T w, T u, T f );
  
  /**
   * \fn vector< T > boundingBox()
   *
   * \brief This method return the bounding box delimiting
   * the region where will be created the shape.
   *
   * \return Vector containing 2 positions with tuple 
   * information the limits of rectangular region ( delta and size ).
   */
  vector< T > boundingBox();
  
  /**
   * \fn Shape< T >* getShape()
   *
   * \brief This method return the shape associated
   *
   * \return The shape level constructed.
   */
  Shape< T >* getShape();

  /**
   * \fn int getTypeShape()
   *
   * \brief This method return the type of shape associated
   *
   * \return Type of shape constructed.
   */
  int getTypeShape();
    
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
  Shape< T >* shape = NULL; ///< Shape of level
  int shapeSaveMode; ///< Block shape 
  vector< T > boundingBoxShape; ///< Bounding box
  int indexLevel; ///< Index of level
  T dimW; ///< Dimension of multiresolution
};

#endif

/**
 * \fn Level( int k, int m, T w, T u, T f, int shapeMode )
 *
 * \brief Default constructor for Level class which will create 
 * a level of Rectangle shape.
 *
 * \param k - Level of fovea
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) to build the fovea
 * \param shapeMode - Feature specification configured (see settings 
 * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
 */
template <typename T>
Level< T >::Level( int k, int m, T w, T u, T f, int shapeMode ){
  shapeSaveMode = shapeMode;
  indexLevel = k;
  dimW = w;
  vector< T > _boundingBox = this->boundingBox( k, m, w, u, f );
  switch ( shapeMode ){
  case _BLOCKS_: {
#ifdef DEBUG
    cout << "Blocks shape actived" << endl;
#endif
    Blocks< T >* block = new Blocks< T >( _boundingBox );
    shape = block;
    break;
  }
  case _RECTANGLE_: {
#ifdef DEBUG
    cout << "Rectangle shape actived" << endl;
#endif
    Rectangle< T >* rectangle = new Rectangle< T >( _boundingBox );
    shape = rectangle;
    break;
  }
  case _POLYGONS_: {
#ifdef DEBUG
    cout << "Polygons shape actived" << endl;
#endif
    Polygons< T >* polygons = new Polygons< T >( _boundingBox, 2 );
    shape = polygons;
      break;
  }
  default:
    cout << "Shape wasn't configured" << endl;
    break;
  }
  boundingBoxShape = _boundingBox;
}


/**
 * \fn ~Level()
 *
 * \brief Default destructor class
 */
template <typename T>
Level< T >::~Level(){
  // It does not need to be implemented
}

/**
 * \fn void updateLevel( int m, T w, T u, T f, int shapeMode )
 *
 * \brief Default constructor for Level class which will create 
 * a level of Rectangle shape.
 *
 * \param m - Number levels of fovea
 * \param w - Size of levels
 * \param u - Size of image
 * \param f - Position (x, y) to build the fovea
 * \param shapeMode - Feature specification configured (see settings 
 * above), where 0 (blocks), 1 (rectangle) or 2 (polygons)
 */
template <typename T>
void
Level< T >::updateLevel( int m, T w, T u, T f, int shapeMode ){
  shapeSaveMode = shapeMode;
  vector< T > _boundingBox = this->boundingBox( indexLevel, m, w, u, f );
  switch ( shapeMode ){
  case _BLOCKS_: {
#ifdef DEBUG
    cout << "Blocks shape actived" << endl;
#endif
    Blocks< T >* block = new Blocks< T >( _boundingBox );
    shape = block;
    break;
  }
  case _RECTANGLE_: {
#ifdef DEBUG
    cout << "Rectangle shape actived" << endl;
#endif
    Rectangle< T >* rectangle = new Rectangle< T >( _boundingBox );
    shape = rectangle;
    break;
  }
  case _POLYGONS_: {
#ifdef DEBUG
    cout << "Polygons shape actived" << endl;
#endif
    Polygons< T >* polygons = new Polygons< T >( _boundingBox, 2 );
    shape = polygons;
    break;
  }
  default:
    cout << "Shape wasn't configured" << endl;
    break;
  }
  boundingBoxShape = _boundingBox;
}

/**
 * \fn void updateLevel( vector< Level< T > > levels  )
 *
 * \brief Default constructor for Level class which will create 
 * a level of Blocks shape.
 *
 * \param levels - List of all previously processed levels
 */
template <typename T>
void
Level< T >::updateLevel( vector< Level< T > > levels ){
  switch ( shapeSaveMode ){
  case _BLOCKS_: {
#ifdef DEBUG
    cout << "Update approach using blocks shape" << endl;
#endif
    vector< Shape< T >* > shapes;
    for ( int l = 0; l < levels.size(); l++ ){
      // Consider only Blocks shape
      shapes.push_back( levels[l].getShape() );
    }
    Blocks< T > s = dynamic_cast<Blocks< T >& >(*shape);
    cout << "Vertices before update" << endl;
    s.printVertices();
    s.breakBlocks( shapes );
    cout << "Vertices after update" << endl;
    s.printVertices();
    shape->updateVertices( s.getVertices() );
    //shape = dynamic_cast< Shape< T >* >(&s);
    break;
  }
  case _RECTANGLE_: {
#ifdef DEBUG
    cout << "Rectangle wasn't configured to update yet" << endl;
#endif
    break;
  }
  case _POLYGONS_: {
#ifdef DEBUG
    cout << "Polygons wasn't configured to update yet" << endl;
#endif
    break;
  }
  default:
    cout << "Shape wasn't configured" << endl;
    break;
  }
}

/**
 * \fn Mat getLevel( Mat img )
 *
 * \brief Method responsable to create an image with 
 * with level dimension
 *
 * \param img - Image will be foveated
 *
 * \return Image that represent the level
 */
template <typename T>
Mat
Level< T >::getLevel( Mat img ){
  // It's necessary cut the image depending of vertices him
  vector< T > boundingBox = shape->getBoundingBox();
  Rect roi((int)boundingBox[0].x, (int)boundingBox[0].y, (int)boundingBox[1].x - 1, (int)boundingBox[1].y - 1);
  Mat croppedImage = img;
  if ( roi.x >= 0 && roi.y >= 0 && roi.width + roi.x < img.cols && roi.height + roi.y < img.rows ) {
    croppedImage = img( roi );
    resize(croppedImage, croppedImage, Size((int)dimW.x, (int)dimW.y), 0, 0, INTER_AREA);
  }
  return croppedImage;
}

/**
 * \fn vector< Mat > getLevelParts( Mat img )
 *
 * \brief Method responsable to create an image with 
 * with parts of level
 *
 * \param img - Image will be foveated
 *
 * \return Image that represent the level
 */
template <typename T>
vector< Mat >
Level< T >::getLevelParts( Mat img ){
  vector< T > boundingBox = shape->getBoundingBox();
  if ( shapeSaveMode == _BLOCKS_ )
    boundingBox = shape->getVertices();
  vector< Mat > output;
  for ( int r = 0; r < boundingBox.size(); r+=2 ){
    Rect roi((int)boundingBox[r].x, (int)boundingBox[r].y, (int)boundingBox[r+1].x - 1, (int)boundingBox[r+1].y - 1);
    if ( roi.x >= 0 && roi.y >= 0 && roi.width + roi.x < img.cols && roi.height + roi.y < img.rows ) {
      Mat croppedImage = img( roi );
      resize(croppedImage, croppedImage, Size((int)dimW.x, (int)dimW.y), 0, 0, INTER_AREA);
      output.push_back( croppedImage );
    }
  }      
  return output;
}

/**
 * \fn vector< T > boundingBox( int k, int m, T w, T u, T f )
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
vector< T > 
Level< T >::boundingBox( int k, int m, T w, T u, T f ){
  T delta = getDelta( k, m, w, u, f ); ///< Delta is the upper left corner of bounding box
  T size = getSize( k, m, w, u ); ///< Size is the dimension between Delta and bottom right corner of bounding box
  vector< T > _boundingBox;  ///< Tuple vector containing delta and size
  _boundingBox.push_back( delta );
  _boundingBox.push_back( size );
  boundingBoxShape = _boundingBox;
  return _boundingBox;
}

/**
 * \fn vector< T > boundingBox()
 *
 * \brief This method return the bounding box delimiting
 * the region where will be created the shape.
 *
 * \return Vector containing 2 positions with tuple 
 * information the limits of rectangular region ( delta and size ).
 */
template <typename T>
vector< T >
Level< T >::boundingBox(){
  return boundingBoxShape; ///< Tuple vector containing delta and size
}

/**
 * \fn Shape< T >* getShape()
 *
 * \brief This method return the shape associated
 *
 * \return The shape level constructed.
 */
template <typename T>
Shape< T >*
Level< T >::getShape(){
  return shape;
}

/**
 * \fn int getTypeShape()
 *
 * \brief This method return the type of shape associated
 *
 * \return Type of shape constructed.
 */
template <typename T>
int
Level< T >::getTypeShape(){
  return shapeSaveMode;
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
  cout << "Delta: ( " << dx << ", " << dy << " ) " << endl;  
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
  int sx = int( ((m * u.x) + (w.x * k) - (k * u.x)) / m );
  int sy = int( ((m * u.y) + (w.y * k) - (k * u.y)) / m );
#ifdef DEBUG
  cout << "Size: ( " << sx << ", " << sy << " ) " << endl;  
#endif
  return T( sx, sy );
}


/** @} */ //end of group class.
