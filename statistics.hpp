/**
 * \file statistics.hpp
 * 
 * \brief This file contains the prototype and implementation of probability
 * calculations and statistics.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at ufrn (dot) edu (dot) br
 *
 * \version 0.1
 * \date January 2019
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

#ifndef STATISTICS_H
#define STATISTICS_H

#include <iostream> //std::cout, std::endl
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector> //std::vector
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/features2d/features2d.hpp"
//#include "fovea.hpp"
#include "multifovea.hpp"
#include "graph.hpp"
#ifdef _OPENMP
#include <omp.h> //#pragma omp parallel for
#endif

/**
 * \defgroup ProjectFovea Project Fovea
 * @{
 */

/**
 * \class Statistics
 *
 * \brief This class implements the Statistics TAD for calculate 
 * statistics of fovea with generic type.
 *
 * \tparam T - Generic representation for type cv::Point
 */
template < typename T > // cv::Point
class Statistics{
public:  
  //
  // Variables
  //

  //
  // Methods
  //
  
  /**
   * \fn void plotProportion( Fovea< T >* fovea, cv::Mat scene, cv::Mat model, int methodDetector,
   * std::vector< cv::KeyPoint > modelKeypoints, cv::Mat modelDescriptors, float threshold ) 
   *
   * \brief Plot the proportion of inliers/(inliers+outliers) in each level.
   *
   * \param fovea - Fovea pointer to be analised
   * \param scene - Image to be foveated
   * \param model - Template/Model to be detected
   * \param method - Feature specification configured (see settings features) 
   * \param modelKeypoints - Model keypoints
   * \param modelDescriptors - Model descriptors
   * \param threshold1, threshold2 - Filtering limits
   */
  void plotProportion( Fovea< T >* fovea, cv::Mat scene, cv::Mat model, int method,
		       std::vector< cv::KeyPoint > modelKeypoints, 
		       cv::Mat modelDescriptors, float threshold1, float threshold2 );
  
  /**
   * \fn double functionFovea( Fovea< T >* fovea, float threshold1, float threshold2 )
   *
   * \brief Calculate the function pondering levels 
   * \f$ f_{nfovea} = \sum_{k = 0}^{k = m + 1} \alpha_{k} Ir_{k} \f$, 
   * where \f$ \alpha \f$ represent the weight of the level and 
   * \f$ Ir \f$ represent the inliers ratio in the level 
   *
   * \param fovea - Fovea pointer to be analised
   * \param threshold1, threshold2 - Filtering limits
   *
   * \return A value that means the fovea's contribution to target detection.
   */
  double functionFovea( Fovea< T >* fovea, float threshold1, float threshold2 );
  
  /**
   * \fn std::vector< double > regionTransformed( T R_t, std::vector< T > R )
   *
   * \brief Calculate the influence of each level. This function is calculated 
   * by following the following idea:
   * \f[
   * weight_{i} = \frac{-R_{i}}{R_{t}} + \sum_{k = 0, k \neq i}^{ m + 1 } \frac{R_{k}}{R_{t}}
   * \f]
   * where \f$ R_{i} \f$ is region size and \f$ R_{t} \f$ is total region.
   *
   * \param R_t - Represent the sum of all regions in space domain
   *        R - Quantity of pixels of transformed region
   *
   * \return returns normalized weights for each level
   */
  std::vector< double > regionTransformed( double R_t, std::vector< double > R );
  
  /**
   * \fn int localGradient( double referencePotential, std::vector< double > potentials )
   *
   * \brief Computes the direction of local gradient
   *
   * \param potentialReference - Indicates the value of the reference potential
   * \param potentials - Contains all the potentials around the reference
   *
   * \return Index related to increased local potential
   */
  int localGradient( double referencePotential, std::vector< double > potentials );

  /**
   * \fn std::vector< T > reduceRegionByLocalGradient( int index, T px, std::vector< T > regionSearched, int configuration )
   *
   * \brief Reduces the region searched according with the highest potential local index
   *
   * \param index - Index related to position of highest potential local
   * \param px - Indicates the position of reference potential analised
   * \param regionSearched - Represent the region with higher probability of find the target
   * \param configuration - Indicates how the potentials were extracted, both in the clockwise direction
   * - If configuration == 0, then potentials extracted from northeast, southeast, southwest and northwest
   * - If configuration == 1, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
   *
   * \return Update region to position fovea
   */
  std::vector< T > reduceRegionByLocalGradient( int index, T px, std::vector< T > regionSearched, int configuration );


  /**
   * \fn T intersectionLocalGradient( std::vector< T > positionPotentialVectorA, 
   *                                  std::vector< T > positionPotentialVectorB )
   *
   * \brief Computes the intersection of local gradients.
   *
   * \param positionPotentialVectorA - Contains the Points (x_ref, y_ref) and (x_max, y_max) for maximum potential from fovea A
   * \param positionPotentialVectorB - Contains the Points (x_ref, y_ref) and (x_max, y_max) for maximum potential from fovea B  
   *
   * \return Point (x, y) that intersect the local gradients.
   */
  T intersectionLocalGradient( std::vector< T > positionPotentialVectorA,
			       std::vector< T > positionPotentialVectorB );

  /**
   * \fn T maximumLikelihoodEstimator( std::vector< T > samples, std::vector< double > potentials, int method )
   *
   * \brief Calculate the Maximum Likelihood Estimator (MLE)
   *
   * \param samples - Contains all points of the features
   * \param potentials - Contains all the potentials of points
   * \param method - Arithmetic mean (0) and weighted average (1)
   *
   * \return Point estimated through MLE
   */
  T maximumLikelihoodEstimator( std::vector< T > samples, std::vector< double > potentials, int method );
  
  /**
   * \fn T trilaterationEstimator( std::vector< T > foveae, std::vector< double > inverseDetectionRate )
   *
   * \brief Calculate the trilateration Estimator
   *
   * \param foveae - Contains all points of the foveae
   *        inverseDetectionRate - Contain inverse detection rate of each fovea
   *
   * \return Point estimated through trilateration
   */
  T trilaterationEstimator( std::vector< T > foveae, std::vector< double > inverseDetectionRate );

  
};

#endif

/**
 * \fn void plotProportion( Fovea< T >* fovea, cv::Mat scene, cv::Mat model, int method,
 * std::vector< cv::KeyPoint > modelKeypoints, cv::Mat modelDescriptors, 
 * float threshold1, float threshold2 ) 
 *
 * \brief Plot the proportion of inliers/(inliers+outliers) in each level.
 *
 * \param fovea - Fovea pointer to be analised
 * \param scene - Image to be foveated
 * \param model - Template/Model to be detected
 * \param methodDetector - Feature specification configured (see settings features) 
 * \param modelKeypoints - Model keypoints
 * \param modelDescriptors - Model descriptors
 * \param threshold1, threshold2 - Filtering limits
 */
template <typename T>
void 
Statistics< T >::plotProportion( Fovea< T >* fovea, cv::Mat scene, cv::Mat model, int method, 
				 std::vector< cv::KeyPoint > modelKeypoints, 
				 cv::Mat modelDescriptors, float threshold1, float threshold2 ){
  std::vector< T > parameters = fovea->getParameters();
  int m = parameters[0].x;
  FILE *file;
  file = fopen ("data/graph3d.dat","w");
  //double maxInliersRatio = 0.0;
  if ( file != NULL ){
    fprintf( file, "#%c #%c #%c \n", 'x', 'y', 'z');
    for ( int i = 0; i < scene.rows; i+=10 ){
      for ( int j = 0; j < scene.cols; j+=10 ){
	fovea->updateFovea( cv::Point( i, j ) );
	fovea->foveatedFeatures( scene, method, MRMF );
	fovea->matching( scene, model, modelKeypoints, modelDescriptors );
	
	fprintf( file, "%d %d ", i, j );
	// Levels
	for ( int k = 0; k < m + 1; k++ ){
	  //if ( k == 0 ) numberMatchesMax = fovea->getNumberMatches( k );
	  if ( ( fovea->getNumberMatches( k ) > threshold1 ) &&
	       ( fovea->getNumberMatches( k ) < threshold2 ) ){
	    fprintf( file, "%f ", fovea->getInliersRatio( k ) );
	    /*if ( maxInliersRatio < fovea->getInliersRatio( k ) ){
	      maxInliersRatio = fovea->getInliersRatio( k );
	      //std::cout << "level = " << k << std::endl; 
	      //std::cout << "( x, y ) = ( " << i << ", " << j << " )" << std::endl;
	    }*/
	  }
	  else{
	    fprintf( file, "%f ", 0.0 );
	  }
	}
	fprintf( file, "\n" );
      }
    }
  }
  fclose( file );
  
  Graph g;
  for ( int k = 0; k < m + 1; k++ ){
    g("set term postscript eps");
    const std::string _file ("data/graph3d");
    g("set output \"data/graph3d_"+std::to_string(k)+".eps\" ");
    // http://lowrank.net/gnuplot/plot3d2-e.html
    g("set view 75, 30, 1, 1");
    g("set xlabel \'x\'");
    g("set ylabel \'y\'");
    //g("set zrange [0:0.2]");
    g("set dgrid3d 30, 30");
    g("set hidden3d");
    g("splot \'"+ _file +".dat\' u 1:2:"+std::to_string(k+3)+" with lines ");
    //g("splot \'"+ _file +".dat\' u 1:2:3 with lines, \'"+ _file +".dat\' u 1:2:4 with lines, \'"+ _file +".dat\' u 1:2:5 with lines ");
  }
}

/**
 * \fn double functionFovea( Fovea< T >* fovea, float threshold1, float threshold2 )
 *
 * \brief Calculate the function pondering levels 
 * \f$ f_{nfovea} = \sum_{k = 0}^{m + 1} \alpha_{k} Ir_{k} \f$, 
 * where \f$ \alpha \f$ represent the weight of the level and 
 * \f$ Ir \f$ represent the inliers ratio in the level 
 *
 * \param fovea - Fovea pointer to be analised
 * \param threshold1, threshold2 - Filtering limits
 *
 * \return A value that means the fovea's contribution to target detection.
 *
 * \note You must filter the ratio of inliers, because few matches can have high inlier rates.
 */
template <typename T>
double 
Statistics< T >::functionFovea( Fovea< T >* fovea, float threshold1, float threshold2 ){
  double function = 0.0;
  std::vector< T > parameters = fovea->getParameters();
  int m = parameters[0].x;
  double totalRegionx = 0.0;
  std::vector< double > regionsx;
  //T totalRegion = T( 0.0, 0.0 );
  //std::vector< T > regions;
  for ( int k = 0; k < m + 1; k++ ){
    if ( ( fovea->getNumberMatches( k ) > threshold1 ) &&
	 ( fovea->getNumberMatches( k ) < threshold2 ) ){
      Level< T > level = fovea->getLevelFromFovea( k );
      std::vector< T > boundingBox = level.boundingBox( k, m, parameters[1], parameters[2], parameters[3] );
      totalRegionx += (1.0/(boundingBox[1].x));
      regionsx.push_back( (1.0/(boundingBox[1].x)) );
      //totalRegion += T( (boundingBox[0].x + boundingBox[1].x), (boundingBox[0].y + boundingBox[1].y) );
      //regions.push_back( T( (boundingBox[0].x + boundingBox[1].x), (boundingBox[0].y + boundingBox[1].y) ) );
    }
  }
  std::vector< double > alpha = regionTransformed( totalRegionx, regionsx );
  //std::vector< double > alpha = regionTransformed( totalRegion, regions );
  for ( int k = 0; k < m + 1; k++ ){
    if ( ( fovea->getNumberMatches( k ) > threshold1 ) &&
	 ( fovea->getNumberMatches( k ) < threshold2 ) ){
      function += alpha[k] * fovea->getInliersRatio( k );
      //std::cout << "level = " << k << " , alfa = " << alpha[k] << " , function = " << fovea->getInliersRatio( k ) << std::endl;
    }
  }
  return function;
}

/**
 * \fn std::vector< double > regionTransformed( T R_t, std::vector< T > R )
 *
 * \brief Calculate the influence of each level. This function is calculated 
 * by following the following idea: 
 * \f$ weight_{i} = \frac{-R_{i}}{R_{t}} + \sum_{k = 0, k \neq i}^{ m + 1 } \frac{R_{k}}{R_{t}} f\$,
 * where \f$ R_{i} \f$ is region size and \f$ R_{t} \f$ is total region.
 *
 * \param R_t - Represent the sum of all regions in space domain
 *        R - Quantity of pixels of transformed region
 *
 * \return returns normalized weights for each level
 */
template <typename T>
std::vector< double >
Statistics< T >::regionTransformed( double R_t, std::vector< double > R ){
  std::vector< double > weights;
  /*for ( int analised = 0; analised < R.size(); analised++ ){
    double complementx = 0.0; double complementy = 0.0;
    for ( int notAnalised = 0; notAnalised < R.size(); notAnalised++ ){
      if ( notAnalised != analised ){
	complementx += ((double)R[notAnalised].x / R_t.x);
	complementy += ((double)R[notAnalised].y / R_t.y);
      }
    }
    complementx -= ((double)R[analised].x / R_t.x);
    complementy -= ((double)R[analised].y / R_t.y);
    weights.push_back( complementx );
  }*/
  
  /*T R_t_complement = T(0, 0);
  for ( int r = 0; r < R.size(); r++ ){
    R_t_complement += R_t - R[r];
  }
  if ( R_t_complement.x != 0 || R_t_complement.y != 0 ){
    for ( int r = 0; r < R.size(); r++ ){
      double wx = double( R_t.x - R[r].x ) / R_t_complement.x;
      double wy = double( R_t.y - R[r].y ) / R_t_complement.y;
      weights.push_back( wx );
    }
  }
  else{
    weights.push_back( 0.0 );
  }*/
  
  for ( int r = 0; r < R.size(); r++ ){
    double wx = double(R[r])/ R_t;
    weights.push_back( wx );
  }
  
  return weights;
}

/**
 * \fn int localGradient( double referencePotential, std::vector< double > potentials )
 *
 * \brief Computes the direction of local gradient
 *
 * \param potentialReference - Indicates the value of the reference potential
 * \param potentials - Contains all the potentials around the reference
 *
 * \return Index related to increased local potential
 */
template <typename T>
int 
Statistics< T >::localGradient( double referencePotential, std::vector< double > potentials ){
  int positionMaxPotential = 0; // Find a position with the higher potential
  for (int i = 1; i < potentials.size(); i++){
    if ( potentials[i] - referencePotential > potentials[positionMaxPotential] - referencePotential )
      positionMaxPotential = i;
  }
  // We can think that's possible to use referencePotential to inform the vector displacement..
  // If referencePotential is low, then the vector displacement is larger
  std::cout << positionMaxPotential;
  return positionMaxPotential;
}

/**
 * \fn std::vector< T > reduceRegionByLocalGradient( int index, T px, std::vector< T > regionSearched, int configuration )
 *
 * \brief Reduces the region searched according with the highest potential local index
 *
 * \param index - Index related to position of highest potential local
 * \param px - Indicates the position of reference potential analised
 * \param regionSearched - Represent the region with higher probability of find the target
 * \param configuration - Indicates how the potentials were extracted, both in the clockwise direction
 * - If configuration == 0, then potentials extracted from northeast, southeast, southwest and northwest
 * - If configuration == 1, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
 *
 * \return Update region to position fovea
 */
template <typename T>
std::vector< T > 
Statistics< T >::reduceRegionByLocalGradient( int index, T px, std::vector< T > regionSearched, int configuration ){
  std::vector< T > updatedRegion = regionSearched;
  T delta = regionSearched[0];
  T size = regionSearched[1];
  switch ( configuration ){
  case 0: // index related to northeast, southeast, southwest and northwest
    if ( index == 0 ){
      updatedRegion[0] = T( delta.x, px.y );
      updatedRegion[1] = T( px.x, size.y - px.y );
    }
    if ( index == 1 ){
      updatedRegion[0] = T( px.x, px.y );
      updatedRegion[1] = T( size.x - px.x, size.y - px.y );
    }
    if ( index == 2 ){
      updatedRegion[0] = T( px.x, delta.y );
      updatedRegion[1] = T( size.x - px.x, px.y );
    }
    if ( index == 3 ){
      //updatedRegion[0] = T( delta.x, delta.y );
      updatedRegion[1] = T( px.x, px.y );
    }
    break;
  case 1: // index related to up (North), right (EAST), down (SOUTH) and left (WEST)
    if ( index == 0 ){
      //updatedRegion[0] = T( delta.x, delta.y ); 
      updatedRegion[1] = T( px.x, size.y );
    }
    if ( index == 1 ){
      updatedRegion[0] = T( delta.x, px.y );
      updatedRegion[1] = T( size.x, size.y - px.y );
    }
    if ( index == 2 ){
      updatedRegion[0] = T( px.x, delta.y );
      updatedRegion[1] = T( size.x - px.x, size.y );
    }
    if ( index == 3 ){
      //updatedRegion[0] = T( delta.x, delta.y );
      updatedRegion[1] = T( size.x, px.y );
    }
    break;
  default:
    std::cout << "This configuration was not configured!" << std::endl;
  }
  return updatedRegion;
}

/**
 * \fn T intersectionLocalGradient( std::vector< T > positionPotentialVectorA, 
 *                                  std::vector< T > positionPotentialVectorB )
 *
 * \brief Computes the intersection of local gradients.
 *
 * \param positionPotentialVectorA - Contains the Points (x_a, y_a) and (x_b, y_b) for maximum potential from fovea A
 *        positionPotentialVectorB - Contains the Points (x_c, y_c) and (x_d, y_d) for maximum potential from fovea B  
 *
 * \return Point (x, y) that intersect the local gradients.
 */
template <typename T>
T
Statistics< T >::intersectionLocalGradient( std::vector< T > positionPotentialVectorA,
					    std::vector< T > positionPotentialVectorB ){
  // Points
  int xa = positionPotentialVectorA[0].x; int ya = positionPotentialVectorA[0].y;
  int xb = positionPotentialVectorA[1].x; int yb = positionPotentialVectorA[1].y;
  int xc = positionPotentialVectorB[0].x; int yc = positionPotentialVectorB[0].y;
  int xd = positionPotentialVectorB[1].x; int yd = positionPotentialVectorB[1].y;
  double x = 0.0; double y = 0.0;
  // angular coefficient
  if ( ( xb - xa != 0 ) || ( xd - xc != 0 ) ){
    double mA = (yb - ya)/(xb - xa);
    double mB = (yd - yc)/(xd - xc);
    std::cout << "mA = " << mA << ", mB = " << mB << std::endl;
    // linear coefficient
    // y = mx - mx_{1} + y_{1} => y = mx + c
    double cA = mA * xa + ya;
    double cB = mB * xc + yc;
    if ( mA == mB ){
      std::cout << "Angular coefficients are equals, then parallel " << std::endl;
      if ( cA != cB )
	std::cout << "and distintes lines" << std::endl;
      if ( cA == cB )
	std::cout << "and coincident lines" << std::endl;
    }
    else{
      if ( mA * mB == -1 ){
	std::cout << "Perpendicular lines" << std::endl;	
	x = ( yc - ya + (mA * xa) - (mB * xc) ) / (mA - mB);
	y = mA * ( x - xa ) + ya;
	
      }
      if ( ( mA == 0 ) || ( mB == 0 ) ){
	std::cout << "Vertical parallel line" << std::endl;
	x = ( yc - ya + (mA * xa) - (mB * xc) ) / (mA - mB);
	y = mA * ( x - xa ) + ya;
      }
      
    }
  }
  else{
    std::cout << "Horizontal parallel lines" << std::endl;
  }
  
  return T(x, y);
}

/**
 * \fn T maximumLikelihoodEstimator( std::vector< T > samples, std::vector< double > potentials, int method )
 *
 * \brief Calculate the Maximum Likelihood Estimator (MLE)
 *
 * \param samples - Contains all points of the features
 * \param potentials - Contains all the potentials of points
 * \param method - Arithmetic mean (0) and weighted average (1)
 *
 * \return Point estimated through MLE
 */
template <typename T>
T
Statistics< T >::maximumLikelihoodEstimator( std::vector< T > samples, std::vector< double > potentials, int method ){
  double x = 0.0; double y = 0.0; double n = 0.0;
  if ( method == 0 ){ // Arithmetic mean
    for (int i = 0; i < samples.size(); i++){
      x += samples[i].x;
      y += samples[i].y;
    }
    x /= samples.size();
    y /= samples.size();
  }
  else{  // weighted average
    for (int i = 0; i < samples.size(); i++){
      x += samples[i].x * potentials[i];
      y += samples[i].y * potentials[i];
      n += potentials[i];
    }
    x /= n;
    y /= n;
  }

  return T(x, y);
}
  
/**
 * \fn T trilaterationEstimator( std::vector< T > foveae, std::vector< double > inverseDetectionRate )
 *
 * \brief Calculate the trilateration Estimator
 *
 * \param foveae - Contains all points of the foveae
 *        inverseDetectionRate - Contain inverse of detection rate for each fovea
 *
 * \return Point estimated through trilateration
 */
template <typename T>
T
Statistics< T >::trilaterationEstimator( std::vector< T > foveae, std::vector< double > inverseDetectionRate ){
  // Points
  int x1 = foveae[0].x; int y1 = foveae[0].y;
  int x2 = foveae[1].x; int y2 = foveae[1].y;
  int x3 = foveae[2].x; int y3 = foveae[2].y;
  // Inverse of detection rate
  double r1 = inverseDetectionRate[0];
  double r2 = inverseDetectionRate[1];
  double r3 = inverseDetectionRate[2];

  double a = (-2 * x1) + (2 * x2);
  double b = (-2 * y1) + (2 * y2);
  double c = (r1*r1) - (r2*r2) - (x1*x1) + (x2*x2) - (y1*y1) + (y2*y2);
  double d = (-2 * x2) + (2 * x3);
  double e = (-2 * y2) + (2 * y3);
  double f = (r2*r2) - (r3*r3) - (x2*x2) + (x3*x3) - (y2*y2) + (y3*y3);

  double den = (a * e) - (b * d);
  double x = ((c * e) - (b * f)) / den;
  double y = ((a * f) - (c * d)) / den;
  
  return T(x, y);
}

/** @} */ //end of group class.
