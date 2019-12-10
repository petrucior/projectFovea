#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "../statistics.hpp"

using namespace cv;

static void on_mouse(int event, int x, int y, int flags, void *_fovea) {
  if ( event == EVENT_MOUSEMOVE ){
    Fovea< Point > *fovea = (Fovea< Point > *) _fovea;
    fovea->updateFovea( Point(x, y ) );
    //std::cout << "( " << x << ", " << y << " )" << std::endl;
  }
}

int main( int argc, char** argv ){
  // Sources
  cv::Mat model = cv::imread(argv[1]);//, CV_LOAD_IMAGE_GRAYSCALE);
  cv::Mat scene = cv::imread(argv[2]);//, CV_LOAD_IMAGE_GRAYSCALE):
  // Multifovea
  cv::Mat multifoveated;
  std::vector< cv::Point > fs;
  // Foveas centralized position
  fs.push_back( cv::Point( 0, 0 ) );
  fs.push_back( cv::Point( 0, 0 ) );
  fs.push_back( cv::Point( 0, 0 ) );
  fs.push_back( cv::Point( 0, 0 ) );
  fs.push_back( cv::Point( 0, 0 ) );
  
  // Foveas simetrically positioned
  /*fs.push_back( cv::Point( -40, 40 ) );
  fs.push_back( cv::Point( 40, 40 ) );
  fs.push_back( cv::Point( 40, -40 ) );
  fs.push_back( cv::Point( -40, -40 ) );*/
  std::vector<cv::Scalar> colors;
  srand (time(NULL)); // Initialize random seed                                                                                                                                                 
  for (int i = 0; i < fs.size(); i++){
    int r = rand() % 256; // 0 - 255                                                                                                                                                            
    int g = rand() % 256; // 0 - 255                                                                                                                                                            
    int b = rand() % 256; // 0 - 255                                                                                                                                                            
    colors.push_back(cv::Scalar(b, g, r));
  }
  
  // Extraction features of model
  cv::Ptr<cv::FeatureDetector> detector;
  cv::Ptr<cv::DescriptorExtractor> descriptor;
  
  // Configuration to _ORB_
  detector = cv::ORB::create();
  descriptor = cv::ORB::create();
  
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 5, cv::Point(120, 120), cv::Point( scene.cols, scene.rows ), fs );
  Multifovea< cv::Point > *fscene = new Multifovea< cv::Point >( 5, cv::Point(120, 120), cv::Point( scene.cols, scene.rows ), fs, REEXECUTION );
  

  // Configuration to _KAZE_
  //detector = cv::KAZE::create();
  //descriptor = cv::KAZE::create();
    
  //Fovea< cv::Point > *fscene = new Fovea< cv::Point >( 4, cv::Point(100, 100), cv::Point( scene.cols, scene.rows ), fs );
  //Multifovea< cv::Point > *fscene = new Multifovea< cv::Point >( 4, cv::Point(60, 60), cv::Point( scene.cols, scene.rows ), fs, REEXECUTION );

  // Extracting model features
  std::vector< cv::KeyPoint > modelKeypoints;
  cv::Mat modelDescriptors;
  detector->detect ( model, modelKeypoints );
  descriptor->compute ( model, modelKeypoints, modelDescriptors );
  
  // Extracting scene features
  multifoveated = fscene->multifoveaLevelsImage( scene, colors );
  fscene->foveatedFeatures( scene, _ORB_, MRMF );
  //int64 t = cv::getTickCount();
  // Matching between model and foveated scene
  fscene->matching( scene, model, modelKeypoints, modelDescriptors );
  //t = cv::getTickCount() - t;
  //std::cout << "time = " << t*1000/cv::getTickFrequency() << " ms ";
  
  Statistics< cv::Point > *s = new Statistics< cv::Point >();
  std::vector< Fovea< cv::Point >* > foveas = fscene->getFoveas();
  for( int f = 0; f < fs.size(); f++ ){
    std::cout << s->functionFovea( foveas[f], 87, 1000 ) << std::endl;
  }
  //imshow("image", multifoveated );
  //waitKey( 0 );

  
  

  /*while ( true ){
    //setMouseCallback( "image", on_mouse, f );
    //multifoveated = m->multifoveatedImage( img );
    multifoveated = m->multifoveaLevelsImage( img, colors );
    imshow("image", multifoveated );
    waitKey( 0 );
    char key = waitKey( 0 );
    if ( key == 'q' ) break;
    if ( key == 'u' )
      m->foveatedFeatures( img, _KAZE_, MRMF );
  }*/
  
  return 0;
}
