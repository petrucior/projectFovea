/**
 * \file statistics.h
 * 
 * \brief This file contains the prototype and implementation of probability
 * calculations and statistics.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date January 2019
 *
 * \copyright
 * Copyright (C) 2016, Petrúcio Ricardo <petrucior@gmail.com>
 * If you use this software for academic purposes, consider citing the related
 * paper: Rafael Beserra Gomes, Bruno Motta de Carvalho, Luiz Marcos Garcia
 * Gonçalves, Visual attention guided features selection with foveated images,
 * Neurocomputing, Volume 120, 23 November 2013, Pages 34-44, ISSN 0925-2312, 
 * http://dx.doi.org/10.1016/j.neucom.2012.10.033.
 *
 * This file is part of foveatedFeatures software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdlib.h>
#include <iostream>
#include <vector>

/**
 * \struct Statistics
 *
 * \brief Struct for calculate statistics.
 */
struct Statistics {

  //
  // Variables
  //

  //
  // Methods
  //

  /**
   * \fn float proportion( int inliers, int outliers )
   *
   * \brief Calculate the proportion of inliers/outliers.
   *
   * \param inliers - Number of inliers
   *        outliers - Number of outliers
   */
  float proportion( int inliers, int outliers );

  /**
   * \fn int factorial( int value )
   *
   * \brief Calculate the function factorial of value.
   *
   * \param value - Number that will be factored
   */
  int factorial( int value );
  
  /**
   * \fn float fieldVision( int m, int k )
   *
   * \brief Calculate the probability of detection by level.
   *
   * \param m - Number of levels minus one
   *        k - Number of current level
   */
  float fieldVision( int m, int k );

   /**
   * \fn float regionTransformed( int R_t, int R_kc )
   *
   * \brief Calculate the influence of each level.
   *
   * \param R_t - Represent the sum of all regions in space domain
   *        R_kc - Complementary quantity of pixels transformed
   */
  float regionTransformed( int R_t, int R_kc );

  /**
   * \fn void localGradient( float referencePotential, std::vector<float> potentials, int configuration )
   *
   * \brief Computes the direction of local gradient.
   *
   * \param potentialReference - Indicates the value of the reference potential
   *        potentials - Contains all the potentials around the reference
   *        configuration - Indicates how the potentials were extracted, both in the clockwise direction
   * - If configuration == 0, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
   * - If configuration == 1, then potentials extracted from northeast, southeast, southwest and northwest
   * - If configuration == 2, then potentials extracted from north, northeast, east, southeast, south, southwest, west and northwest
   */
  void localGradient( float referencePotential, std::vector<float> potentials, int configuration );

  /**
   * \fn cv::Point intersectionLocalGradient( std::vector< cv::Point > positionPotentialVectorA, 
   *                                     std::vector< cv::Point > positionPotentialVectorB )
   *
   * \brief Computes the intersection of local gradients.
   *
   * \param positionPotentialVectorA - Contains the Points (x_ref, y_ref) and (x_max, y_max) for maximum potential from fovea A
   *        positionPotentialVectorB - Contains the Points (x_ref, y_ref) and (x_max, y_max) for maximum potential from fovea B  
   *
   * \return Point (x, y) that intersect the local gradients.
   */
  cv::Point intersectionLocalGradient( std::vector< cv::Point > positionPotentialVectorA,
				       std::vector< cv::Point > positionPotentialVectorB );

  /**
   * \fn cv::Point maximumLikelihoodEstimator( std::vector< cv::Point > samples )
   *
   * \brief Calculate the Maximum Likelihood Estimator (MLE)
   *
   * \param samples - Contains all points of the features
   *
   * \return Point estimated through MLE
   */
  cv::Point maximumLikelihoodEstimator( std::vector< cv::Point > samples );

  /**
   * \fn cv::Point trilaterationEstimator( std::vector< cv::Point > foveae, std::vector< float > inverseDetectionRate )
   *
   * \brief Calculate the trilateration Estimator
   *
   * \param foveae - Contains all points of the foveae
   *        inverseDetectionRate - Contain inverse detection rate of each fovea
   *
   * \return Point estimated through trilateration
   */
  cv::Point trilaterationEstimator( std::vector< cv::Point > foveae, std::vector< float > inverseDetectionRate );

  
};

#endif


/**
 * \fn float proportion( int inliers, int outliers )
 *
 * \brief Calculate the proportion of inliers/outliers.
 *
 * \param inliers - Number of inliers
 *        outliers - Number of outliers
 */
float
Statistics::proportion( int inliers, int outliers ){
  return inliers/(inliers + outliers);
}


/**
 * \fn int factorial( int value )
 *
 * \brief Calculate the function factorial of value.
 *
 * \param value - Number that will be factored.
 */
int
Statistics::factorial( int value ){
  return ( (value == 0) || (value == 1) ) ? 1 : factorial(value - 1) * value;
}

/**
 * \fn float fieldVision( int m, int k )
 *
 * \brief Calculate influence of each level.
 *
 * \param m - Number of levels minus one
 *        k - Number of current level
 */
float
Statistics::fieldVision( int m, int k ){
  return ( (m + 1) - k )/( factorial(m + 1) );
}


/**
 * \fn float regionTransformed( int R_t, int R_kc )
 *
 * \brief Calculate the influence of each level.
 *
 * \param R_t - Represent the sum of all regions in space domain
 *        R_kc - Complementary quantity of pixels transformed
 */
float
Statistics::regionTransformed( int R_t, int R_kc ){
  return ( R_kc / R_t );
}
  
/**
 * \fn void localGradient( float referencePotential, std::vector<float> potentials, int configuration )
 *
 * \brief Computes the direction of local gradient.
 *
 * \param potentialReference - Indicates the value of the reference potential
 *        potentials - Contains all the potentials around the reference
 *        configuration - Indicates how the potentials were extracted, both in the clockwise direction
 * - If configuration == 0, then potentials extracted from up (North), right (EAST), down (SOUTH) and left (WEST)
 * - If configuration == 1, then potentials extracted from northeast, southeast, southwest and northwest
 * - If configuration == 2, then potentials extracted from north, northeast, east, southeast, south, southwest, west and northwest
 */
void
Statistics::localGradient( float referencePotential, std::vector<float> potentials, int configuration ){
  int positionMaxPotential = 0;
  for (int i = 0; i < potentials.size(); i++){
    if ( potencials[i] - referecePotential > potencials[positionMaxPotential] - referecePotential )
      positionMaxPotential = i;
  }
  // We can think that's possible to use referencePotential to inform the vector displacement..
  // If referencePotential is low, then the vector displacement is larger
  std::cout << positionMaxPotential;
}

/**
 * \fn cv::Point intersectionLocalGradient( std::vector< cv::Point > positionPotentialVectorA, 
 *                                     std::vector< cv::Point > positionPotentialVectorB )
 *
 * \brief Computes the intersection of local gradients.
 *
 * \param positionPotentialVectorA - Contains the Points (x_a, y_a) and (x_b, y_b) for maximum potential from fovea A
 *        positionPotentialVectorB - Contains the Points (x_c, y_c) and (x_d, y_d) for maximum potential from fovea B  
 *
 * \return Point (x, y) that intersect the local gradients.
 */
cv::Point
Statistics::intersectionLocalGradient( std::vector< cv::Point > positionPotentialVectorA,
				       std::vector< cv::Point > positionPotentialVectorB ){
  // Points
  int xa = positionPotentialVectorA[0].x; int ya = positionPotentialVectorA[0].y;
  int xb = positionPotentialVectorA[1].x; int yb = positionPotentialVectorA[1].y;
  int xc = positionPotentialVectorB[0].x; int yc = positionPotentialVectorB[0].y;
  int xd = positionPotentialVectorB[1].x; int yd = positionPotentialVectorB[1].y;
  float x = ( (xb - ya)*((-(yd - yc)*xc) + ((xd - yc)*yc)) - (xd - yc)*((-(yb - ya)*xa) + ((xb - ya)*ya)) ) / ( ((xd - xc)*(yb - ya)) - ((xb - ya)*(yd - yc)) );
  float y = ( ( ((yd - yc)*x) - ((yd - yc)*xc) + ((xd - yc)*yc) )/(xd - yc) );
  return cv::Point(x, y);
}

/**
 * \fn cv::Point maximumLikelihoodEstimator( std::vector< cv::Point > samples )
 *
 * \brief Calculate the Maximum Likelihood Estimator (MLE)
 *
 * \param samples - Contains all points of the features
 *
 * \return Point estimated through MLE
 */
cv::Point
Statistics::maximumLikelihoodEstimator( std::vector< cv::Point > samples ){
  float x = 0.0; float y = 0.0;
  for (int i = 0; i < samples.size(); i++){
    x += samples[i].x;
    y += samples[i].y;
  }
  x /= samples.size();
  y /= samples.size();
  return cv::Point(x, y);
}
  
/**
 * \fn cv::Point trilaterationEstimator( std::vector< cv::Point > foveae, std::vector< float > inverseDetectionRate )
 *
 * \brief Calculate the trilateration Estimator
 *
 * \param foveae - Contains all points of the foveae
 *        inverseDetectionRate - Contain inverse of detection rate for each fovea
 *
 * \return Point estimated through trilateration
 */
cv::Point
Statistics::trilaterationEstimator( std::vector< cv::Point > foveae, std::vector< float > inverseDetectionRate ){
  // Points
  int x1 = foveae[0].x; int y1 = foveae[0].y;
  int x2 = foveae[1].x; int y2 = foveae[1].y;
  int x3 = foveae[2].x; int y3 = foveae[2].y;
  // Inverse of detection rate
  int r1 = inverseDetectionRate[0];
  int r2 = inverseDetectionRate[1];
  int r3 = inverseDetectionRate[2];

  int a = (-2 * x1) + (2 * x2);
  int b = (-2 * y1) + (2 * y2);
  int c = (r1*r1) - (r2*r2) - (x1*x1) + (x2*x2) - (y1*y1) + (y2*y2);
  int d = (-2 * x2) + (2 * x3);
  int e = (-2 * y2) + (2 * y3);
  int f = (r2*r2) - (r3*r3) - (x2*x2) + (x3*x3) - (y2*y2) + (y3*y3);

  int d = (a * e) - (b * d);
  float x = ((c * e) - (b * f)) / d;
  float y = ((a * f) - (c * d)) / d;

  return cv::Point(x, y);
}
