/**
 * \file detections.h
 * 
 * \brief This file contains the prototype and implementation of strategies
 * to move the foveae.
 *
 * \author
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n 
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date February 2019
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

#ifndef DETECTIONS_H
#define DETECTIONS_H

#include <stdlib.h>
#include <iostream>
#include <vector>
#include "statistics.h"

/**
 * \struct Detections
 *
 * \brief Struct of strategies to control the foveae.
 */
struct Detections {
  
  //
  // Variables
  //

  //
  // Methods
  //

  /**
   * \fn cv::Point lastPosition( cv::Point lastPosition )
   *
   * \brief Detection strategy that holds the last position. 
   *
   * \param lastPosition - Last position used to detect
   *
   * \return Point estimated through last position
   */
  cv::Point lastPosition( cv::Point lastPosition );

  /**
   * \fn cv::Point nFrames( std::vector< cv::Point > samples )
   *
   * \brief Detection strategy that holds the n last frames using MLE, but 
   * is possible to use kalman filter. 
   *
   * \param samples - vector with position of all keypoints of n last frames
   *
   * \return Point estimated through n last position
   */
  cv::Point nFrames( std::vector< cv::Point > samples );

  /**
   * \fn void disableFovea( bool &fovea )
   *
   * \brief Detection strategy that disable the fovea. 
   *
   * \param fovea - Pointer to enable/disable fovea
   */
  void disableFovea( bool &fovea );

  /**
   * \fn void increaseGrowthFactor( int &wx, int &wy, float multiplier )
   *
   * \brief Detection strategy that increase growth factor. 
   *
   * \param wx - Pointer to increase x axis the fovea
   *        wy - Pointer to increase y axis the fovea
   * multiplier - Multiplier that will increase the growth factor
   */
  void increaseGrowthFactor( int &wx, int &wy, float multiplier );

  /**
   * \fn cv::Point saliencyMap( cv::Mat img )
   *
   * \brief Detection strategy that use saliency map to control
   * foveae. 
   *
   * \param img - Source image
   *
   * \return Point estimated through saliency map
   */
  cv::Point saliencyMap( cv::Mat img );
  
};

#endif

/**
 * \fn cv::Point lastPosition( cv::Point lastPosition )
 *
 * \brief Detection strategy that holds the last position. 
 *
 * \param lastPosition - Last position used to detect
 *
 * \return Point estimated through last position
 */
cv::Point
Detections::lastPosition( cv::Point lastPosition ){
  return lastPosition;
}

/**
 * \fn cv::Point nFrames( std::vector< cv::Point > samples )
 *
 * \brief Detection strategy that holds the n last frames using MLE, but 
 * is possible to use kalman filter.
 *
 * \param samples - vector with position of all keypoints of n last frames
 *
 * \return Point estimated through n last position
 */
cv::Point
Detections::nFrames( std::vector< cv::Point > samples ){
  Statistics statistic = new Statistics();
  return statistic.maximumLikelihoodEstimator( samples );
}

/**
 * \fn void disableFovea( bool &fovea )
 *
 * \brief Detection strategy that disable the fovea. 
 *
 * \param fovea - Pointer to enable/disable fovea
 * If fovea == 0, then the fovea is disable
 */
void
Detections::disableFovea( bool &fovea ){
  fovea = 0; // Disabling the fovea
}

/**
 * \fn void increaseGrowthFactor( int &wx, int &wy, float multiplier )
 *
 * \brief Detection strategy that increase growth factor. 
 *
 * \param wx - Pointer to increase x axis the fovea
 *        wy - Pointer to increase y axis the fovea
 * multiplier - Multiplier that will increase the growth factor
 */
void 
Detections::increaseGrowthFactor( int &wx, int &wy, float multiplier ){
  wx *= multiplier;
  wy *= multiplier;
}

/**
 * \fn cv::Point saliencyMap( cv::Mat img )
 *
 * \brief Detection strategy that use saliency map to control
 * foveae. 
 *
 * \param img - Source image
 *
 * \return Point estimated through saliency map
 */
cv::Point
Detections::saliencyMap( cv::Mat img ){
  cv::Mat saliencyMap;
  Ptr<Saliency> saliencyAlgorithm = Saliency::create("SPECTRAL_RESIDUAL");
  saliencyAlgorithm.computeSaliency(img, saliencyMap);
  // SaliencyMap contém o mapa de saliência da imagem passada como argumento
  // a questão é detectar em qual parte do mapa de saliência se encontra o
  // ponto de maior saliência da imagem
  return cv::Point(0, 0);
}


