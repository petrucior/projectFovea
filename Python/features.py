#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file feature.py

\brief This file contains the prototype of feature.

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date February 2020
 
This file is part of projectFovea software.
This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 2 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import cv2
from parameters import Parameters
from level import Level

class Feature :
    '''
    \class Feature
    
    \brief This class implements the Feature TAD
    '''
    
    def __init__( self, levels, parameters ):
        '''
        \fn Feature( levels, parameters )
        
        \brief Constructor default
        
        \param levels - set of levels to build the fovea
        \param parameters - Parameters of fovea structure
        '''
        assert( len(levels) > 0 ) # Feature extraction only from a fovea structure
        # SURF parameters
        hessianThreshold = parameters.hessianThreshold
        nOctaves = parameters.nOctaves
        nOctaveLayers = 3
        extended = false
        upright = false
        # Creating SURF object
        surf = cv2..xfeatures2d_SURF( hessianThreshold, nOctaves, nOctaveLayers, extended, upright )
        self.kp = [] # cleaning keypoints
        self.des = [] # cleaning descriptors
        for k in range( 0, parameters.m + 1 ):
            kpAux, desAux = surf.detectAndCompute( levels[k], None )
            self.kp.append( kpAux )
            self.des.append( desAux )

    
    def showFeatureByLevel( self, levels, parameters ):
        '''
        \fn show( levels, parameters )
        
        \brief Display the multiresolution structure
        
        \param levels - set of levels to build the fovea
        \param parameters - Parameters of fovea structure
        '''
        assert( len(levels) > 0 ) # Feature extraction only from a fovea structure
        for k in range( 0, parameters.m + 1 ):
            imgLevelKp = cv2.drawKeypoints( levels[k], self.kp[k], None, (255,0,0), 4 ) 
            imshow( "level "+ k +" ", imgLevelKp )
            waitKey( 0 )
        
    
#How to instantiate and use this class
#if __name__ == '__main__':
#    params = Parameters('params.yaml')
#    img = cv2.imread('../../box_in_scene.png')
#    rows,cols = img.shape[:2]
#    u = [ cols, rows ]
#    fovea = Fovea( u, params )
#    color = (255, 255, 255) # B G R
#    output = fovea.foveatedImage( img, params, color )
#    cv2.imshow( "foveated image", output )
#    cv2.waitKey( 0 )
