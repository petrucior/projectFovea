/*#include <iostream>
#include <string>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../statistics.hpp"*/

#include <iostream>
#include <stdlib.h>
#include <string>
#include "opencv2/core.hpp"
#include "opencv2/calib3d.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/xfeatures2d.hpp"
#include "../statistics.hpp"

using namespace cv;

static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  if ( event == EVENT_MOUSEMOVE ){
    Fovea< Point2f > *fovea = (Fovea< Point2f > *) _fovea;
    fovea->updateFovea( Point2f(x, y ) );
    std::cout << "( " << x << ", " << y << " )" << std::endl;
  }
}

int main( int argc, char** argv ){
  
  // Testing constructor and print functions
  std::cout << argv[1] << std::endl;
  cv::Mat model = cv::imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  cv::Mat scene = cv::imread(argv[2]);//, CV_LOAD_IMAGE_GRAYSCALE):
  cv::Mat sceneFoveated;
  Mat multifoveated = scene;
  std::vector< cv::Point2f > fs, fr;
  /*int quantFoveas = 10;
  // Foveas randomly positioned
  for( int qf = 0; qf < quantFoveas; qf++ ){
    float x = rand() % scene.cols/2;
    float y = rand() % scene.rows/2;
    if ( x < scene.cols/2 ) { x = (scene.cols/2 - x); } else{ x = (x - scene.cols/2); }
    if ( y < scene.rows/2 ) { y = (scene.rows/2 - y); } else{ y = (y - scene.rows/2); }
    fs.push_back( cv::Point2f(x, y) );
  }*/
  // Foveas simetrically positioned
  fs.push_back( cv::Point2f( -150, 70.0 ) );
  fs.push_back( cv::Point2f( -80, 70.0 ) );
  fs.push_back( cv::Point2f( -115, -20.0 ) );
  fs.push_back( cv::Point2f( -40.0, 40.0 ) );
  fs.push_back( cv::Point2f( 40.0, 40.0 ) );
  fs.push_back( cv::Point2f( 40.0, -40.0 ) );
  fs.push_back( cv::Point2f( -40.0, -40.0 ) );
  fs.push_back( cv::Point2f( -80.0, 80.0 ) );
  fs.push_back( cv::Point2f( 80.0, 80.0 ) );
  fs.push_back( cv::Point2f( 80.0, -80.0 ) );
  fs.push_back( cv::Point2f( -80.0, -80.0 ) );
  fs.push_back( cv::Point2f( -90.0, 50.0 ) );
  fs.push_back( cv::Point2f( 90.0, 50.0 ) );
  fs.push_back( cv::Point2f( 90.0, -50.0 ) );
  fs.push_back( cv::Point2f( -90.0, -50.0 ) );
  std::vector< cv::Scalar > colors;
  srand (time(NULL)); // Initialize random seed                                                                                                                                                 
  for (int i = 0; i < fs.size(); i++){
    int r = rand() % 256; // 0 - 255                                                                                                                                                            
    int g = rand() % 256; // 0 - 255                                                                                                                                                            
    int b = rand() % 256; // 0 - 255                                                                                                                                                            
    colors.push_back( cv::Scalar(b, g, r) );
  }
  // Extraction features of model
  cv::Ptr<cv::FeatureDetector> detector;
  cv::Ptr<cv::DescriptorExtractor> descriptor;

  // -----------------------
  // Configuration to _SURF_
  // -----------------------
  detector = cv::xfeatures2d::SURF::create();
  descriptor = cv::xfeatures2d::SURF::create();
  Multifovea< cv::Point2f > *m = new Multifovea< cv::Point2f >( 5, cv::Point2f(80, 80), cv::Point2f( scene.cols, scene.rows ), fs, _REEXECUTION_, _BLOCKS_ );
    
  // -----------------------
  // Configuration to _ORB_
  // -----------------------
  //detector = cv::ORB::create();
  //descriptor = cv::ORB::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(120, 120), cv::Point( scene.cols, scene.rows ), fs );
  
  // -----------------------
  // Configuration to _KAZE_
  // -----------------------
  //detector = cv::KAZE::create();
  //descriptor = cv::KAZE::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 4, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
    
  // -----------------------
  // Configuration to _AKAZE_
  // -----------------------
  //detector = cv::AKAZE::create();
  //descriptor = cv::AKAZE::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
  
  // -----------------------
  // Configuration to _BRISK_
  // -----------------------
  //detector = cv::BRISK::create();
  //descriptor = cv::BRISK::create();
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
    
  
  // -----------------------
  // Configuration to Model
  // -----------------------
  std::vector< cv::KeyPoint > modelKeypoints;
  cv::Mat modelDescriptors;
  detector->detect ( model, modelKeypoints );
  descriptor->compute ( model, modelKeypoints, modelDescriptors );
  
  Statistics< cv::Point2f > *s = new Statistics< cv::Point2f >();
  //Statistics< cv::Point > *s = new Statistics< cv::Point >();
  //s->plotProportion( fscene, scene, model, _SURF_, modelKeypoints, modelDescriptors, 0.0, 1000.0 );
  
  cv::namedWindow("sceneFoveated", 1);
  while ( true ){
    /*setMouseCallback( "sceneFoveated", on_mouse, fscene );
    sceneFoveated = fscene->foveatedImage( scene, cv::Scalar(255, 255, 255) );
    imshow("sceneFoveated", sceneFoveated );
    drawKeypoints( model, modelKeypoints, model, Scalar::all(-1), DrawMatchesFlags::DEFAULT );
    imshow("modelFoveated", model );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'e' ){
      fscene->foveatedFeatures( scene, _SURF_, MRMF );
      fscene->matching( scene, model, modelKeypoints, modelDescriptors );
      
      std::cout << s->functionFovea( fscene, 15.0, 1000.0 ) << std::endl;
      
    }*/
    
    //setMouseCallback( "sceneFoveated", on_mouse, (m->getFoveas())[0]);
    multifoveated = m->multifoveatedImage( scene );
    //multifoveated = m->multifoveaLevelsImage( scene, colors );
    imshow("sceneFoveated", multifoveated );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'e' ){
      m->foveatedFeatures( scene, _SURF_, MRMF );
      m->matching( scene, model, modelKeypoints, modelDescriptors );
      std::vector< double > potentials = s->functionMultiFovea( m, 15.0, 1000.0 );
      std::vector< cv::Point2f > points;
      std::vector< double > invpotentials;
      std::vector< cv::Point2f > pointsMulti;
      std::vector< double > invpotentialsMulti;
      for ( int f = 0; f < fs.size(); f++ ){
	std::cout << "potentials: " << potentials[f];
	fr.push_back( cv::Point2f( scene.cols/2 + fs[f].x, scene.rows/2 + fs[f].y ) );
	if ( f < 3 ){
	  invpotentials.push_back( ((1.0 - potentials[f]) * 0.001 )/ 0.0001 );
	  points.push_back( cv::Point2f( scene.cols/2 + fs[f].x, scene.rows/2 + fs[f].y ) );
	  std::cout << ", invpotentials: " << invpotentials[f];
	}
	if ( f < 8 ){
	  invpotentialsMulti.push_back( ((1.0 - potentials[f]) * 0.001 )/ 0.0001 );
	  pointsMulti.push_back( cv::Point2f( scene.cols/2 + fs[f].x, scene.rows/2 + fs[f].y ) );
	}
	
	std::cout << std::endl;
      }
      
      cv::Point2f pointEstimated;
      pointEstimated = s->maximumLikelihoodEstimator( fr, potentials, 1 );
      std::cout << "Maximum Likelihood Estimator ( MLE )" << std::endl;
      std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;
      cv::circle( multifoveated, cv::Point( (int)pointEstimated.x /*+ scene.cols/2*/, (int)pointEstimated.y /*+ scene.rows/2*/  ), 4, cv::Scalar(0, 0, 255 ), -1, 8, 0);

      pointEstimated = s->trilaterationEstimator(points, invpotentials);
      std::cout << "Trilateration" << std::endl;
      std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;
      cv::circle( multifoveated, cv::Point( (int)pointEstimated.x /*+ scene.cols/2*/, (int)pointEstimated.y /*+ scene.rows/2*/  ), 4, cv::Scalar(0, 255, 255 ), -1, 8, 0);

      pointEstimated = s->baricentricCoordinates( fr, potentials );
      std::cout << "Baricentric Coordinates" << std::endl;
      std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;

      pointEstimated = s->multilateration(pointsMulti, invpotentialsMulti);
      std::cout << "Multiteration" << std::endl;
      std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;
      //cv::circle( multifoveated, cv::Point( (int)pointEstimated.x /*+ scene.cols/2*/, (int)pointEstimated.y /*+ scene.rows/2*/  ), 4, cv::Scalar(0, 255, 255 ), -1, 8, 0);

      /*points.clear();
      invpotentials.clear();
      points.push_back( cv::Point2f( 4, 8 ) );
      points.push_back( cv::Point2f( 9, 6 ) );
      points.push_back( cv::Point2f( 5, 5 ) );
      invpotentials.push_back( 2.82 );
      invpotentials.push_back( 3.0 );
      invpotentials.push_back( 1.41 );
      pointEstimated = s->multilateration( points, invpotentials );
      std::cout << "Multilateration" << std::endl;
      std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;*/

      /*
      vector< vector< double > > a ( 3, vector< double >(3) );
      a[0][0] = 2; a[0][1] = 3; a[0][2] = 4;
      a[1][0] = 7; a[1][1] = 1; a[1][2] = -3;
      a[2][0] = 4; a[2][1] = 3; a[2][2] = -2;
      vector< vector< double > > b ( 3, vector< double >(1) );
      b[0][0] = -1;
      b[1][0] = 5;
      b[2][0] = 3;
      vector< vector< double > > ab = s->multiplication( a, b );
      vector< vector< double > > at = s->transposed( a );
      vector< vector< double > > w ( 3, vector< double >(3) );
      w[0][0] = 13; w[0][1] = 24; w[0][2] = 11;
      w[1][0] = 22; w[1][1] = 47; w[1][2] = 33;
      w[2][0] = 42; w[2][1] = 38; w[2][2] = -6;
      vector< vector< double > > ainv = s->inverse( w );
      for ( int i = 0; i < ainv.size(); i++ ){
	for ( int j = 0; j < ainv[0].size(); j++ )
	  std::cout << ainv[i][j] << "    ";
	std::cout << endl;
      }
      */
    }
    
  }
  
  /*double referencePotential = 0.8;
  std::vector< double > potentials;  
  potentials.push_back( 0.5 );
  potentials.push_back( 0.7 );
  potentials.push_back( 0.9 );
  potentials.push_back( 0.3 );
  int index = s->localGradient( referencePotential, potentials );
  cv::Point px = cv::Point( 10, 10 );
  std::vector< cv::Point2f > regionSearched, updatedRegion;
  regionSearched.push_back( cv::Point2f( 0, 0 ) );
  regionSearched.push_back( cv::Point2f( 512, 512 ) );
  int configuration = 0;
  updatedRegion = s->reduceRegionByLocalGradient( index, px, regionSearched, configuration );
  std::cout << std::endl;
  for( int i = 0; i < regionSearched.size(); i++ )
    std::cout << "( " << updatedRegion[i].x << ", " << updatedRegion[i].y << " ) " << std::endl;
  */

  /*std::vector< cv::Point2f > pointA, pointB;
  pointA.push_back( cv::Point2f( 1, 0 ) );
  //pointA.push_back( cv::Point2f( -1, -1 ) );
  pointA.push_back( cv::Point2f( 0, 0 ) );
  //pointB.push_back( cv::Point2f( 10, 0 ) );
  //pointB.push_back( cv::Point2f( 8, 2 ) );
  pointB.push_back( cv::Point2f( 9, 1 ) );
  pointB.push_back( cv::Point2f( 11, -1 ) );
  cv::Point2f pointEstimated = s->intersectionLocalGradient(pointA, pointB);
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;
  */
  
  /*std::vector< cv::Point2f > points;
  points.push_back( cv::Point2f( 4, 8 ) );
  points.push_back( cv::Point2f( 9, 6 ) );
  points.push_back( cv::Point2f( 5, 5 ) );
  points.push_back( cv::Point2f( 0, 0 ) );
  points.push_back( cv::Point2f( 7, 7 ) );
  // point to be estimated is (6, 6)
  std::vector< double > potentials;
  
  //      pos ----- potential  \  complement ----- radius
  // ref 6, 6 -----    1       \   0.0001    ----- 0.001
  //     3, 3 -----    0.5     \   0.5       ----- 5
  //     4, 8 -----    0.718   \   0.282     ----- 2.82
  //     9, 6 -----    0.7     \   0.3       ----- 3
  //     5, 5 -----    0.859   \   0.141     ----- 1.41
  //     0, 0 -----    0       \   1         ----- 8.48
  //     7, 7 -----    0.859   \   0.141     ----- 1.41
  
  potentials.push_back(0.718);
  potentials.push_back(0.7);
  potentials.push_back(0.859);
  potentials.push_back(0);
  potentials.push_back(0.859);
  cv::Point2f pointEstimated = s->maximumLikelihoodEstimator( points, potentials, 1 );
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;*/
  
  /*std::vector< cv::Point2f > points;
  points.push_back( cv::Point2f( 4, 8 ) );
  points.push_back( cv::Point2f( 9, 6 ) );
  points.push_back( cv::Point2f( 5, 5 ) );
  //points.push_back( cv::Point2f( 0, 0 ) );
  //points.push_back( cv::Point2f( 7, 7 ) );
  // point to be estimated is (6, 6)
  std::vector< double > invpotentials;
  
  //      pos ----- potential  \  complement ----- radius
  // ref 6, 6 -----    1       \   0.0001    ----- 0.001
  //     3, 3 -----    0.5     \   0.5       ----- 5
  //     4, 8 -----    0.718   \   0.282     ----- 2.82
  //     9, 6 -----    0.7     \   0.3       ----- 3
  //     5, 5 -----    0.859   \   0.141     ----- 1.41
  //     0, 0 -----    0       \   1         ----- 8.48
  //     7, 7 -----    0.859   \   0.141     ----- 1.41
  
  invpotentials.push_back(2.82);
  invpotentials.push_back(3);
  invpotentials.push_back(1.41);
  //invpotentials.push_back(8.48);
  //invpotentials.push_back(1.41);
  cv::Point2f pointEstimated = s->trilaterationEstimator(points, invpotentials);
  std::cout << "( " << pointEstimated.x << ", " << pointEstimated.y << " ) " << std::endl;*/

  
  
  
  return 0;
  
}
