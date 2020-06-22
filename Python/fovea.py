#!/usr/bin/pythonA
# -*- coding: utf-8 -*-

'''
\file fovea.py

\brief This file contains the prototype of fovea.

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date January 2020
 
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
from features import Features

class Fovea :
    '''
    \class Fovea
    
    \brief This class implements the Fovea TAD
    '''
    
    def __init__( self, u, parameters ):
        '''
        \fn Fovea( u, parameters )
        
        \brief Constructor default
    
        \param u - Size of image
        \param parameters - Parameters of fovea structure
        '''
        assert( len(parameters.f) == 2 )
        
        # Updating the image size parameter ( U )
        parameters.updateParameterSizeImage( u )

        self.levels = [] # cleaning levels
        for k in range( 0, parameters.m + 1 ):
            level = Level( k, parameters )
            self.levels.append( level )

    def updateFovea( self, parameters ):
        '''
        \fn updateFovea( parameters )

        \brief This method update the fovea structure

        \param parameters - Parameters of fovea structure
        '''
        self.levels = [] # cleaning levels
        for k in range( 0, parameters.m + 1 ):
            level = Level( k, parameters )
            self.levels.append( level )

    def foveatedFeatures( self, img, parameters ):
        '''
        \fn foveatedFeatures( parameters )

        \brief Compute keypoints and descriptors of levels

        \param img - Image to be foveated
        \param parameters - Parameters of fovea structure
        '''
        if ( parameters.typeFovea == 0 ): # MRMF
            self.features = Features( img, self.levels, parameters )
        else: # MMF
            print("coming soon . . .")


    def foveatedMatching( self, model, parameters ):
        '''
        \fn foveatedMatching( model, parameters )

        \brief Maching between multiresolution domain and model 

        \param model - Template that represent the visual stimulli
        \param parameters - Parameters of fovea structure
        '''
        self.features.matchingFeatures( model, parameters )
    
    
    def foveatedImage( self, img, parameters ):
        '''
        \fn foveatedImage( img, parameters )

        \brief This function builds the focused image.

        \param img - Image to be foveated
        \param parameters - Parameters of fovea structure
        '''
        p = parameters
        color = ( p.colors[0], p.colors[1], p.colors[2] )
        imgFoveated = img.copy()
        for k in range( 0, p.m + 1 ):
            imgLevel = self.levels[k].getLevel( img, p )
            pi = self.mapLevel2Image( k, p, [0, 0] )
            pf = self.mapLevel2Image( k, p, p.w )
            if ( k < p.m ):
                imgFoveated[ int(pi[1]):int(pf[1]), int(pi[0]):int(pf[0]) ] = cv2.resize( imgLevel, ( int(pf[0]) - int(pi[0]), int(pf[1]) - int(pi[1]) ) )
            else:
                imgFoveated[ int(pi[1]):int(pf[1]), int(pi[0]):int(pf[0]) ] = imgLevel
            imgFoveated = cv2.rectangle(imgFoveated, (int(pi[0]), int(pi[1])), (int(pf[0]-1), int(pf[1]-1)), color )
        return imgFoveated
    
    
    def mapLevel2Image( self, k, parameters, px ):
        '''
        \fn mapLevel2Image( k, parameters )

        \brief Calculates the position of pixel on the level to image.

        \param k - Level of fovea
        \param parameters - Parameters of fovea structure
        \param px - Pixel (x, y) that we want to map.
        '''
        p = parameters
        pointx = ( (k * p.w[0]) * (p.u[0] - p.w[0]) + (2 * k * p.w[0] * p.f[0]) + (2 * px[0]) * ( (p.m * p.u[0]) - (k * p.u[0]) + (k * p.w[0]) ) )/ (2 * p.m * p.w[0]);
        pointy = ( (k * p.w[1]) * (p.u[1] - p.w[1]) + (2 * k * p.w[1] * p.f[1]) + (2 * px[1]) * ( (p.m * p.u[1]) - (k * p.u[1]) + (k * p.w[1]) ) )/ (2 * p.m * p.w[1]);
        return [ pointx, pointy ]

    
    def saveFoveatedImage( self, name, img, parameters ):
        '''
        \fn saveFoveatedImage()

        \brief Saves all levels to a PNG file

        \param name - Image name
        \param img - Image file
        \param parameters - Parameters of fovea structure
        '''
        print('saving foveated image')
        cv2.imwrite( "midia/" + name + ".png", img )

    
    def saveLevels( self, name, img, parameters ):
        '''
        \fn saveLevels()

        \brief Saves all levels to a PNG file

        \param name - Image name
        \param img - Image file
        \param parameters - Parameters of fovea structure
        '''
        print('saving images')
        for k in range(0, len(self.levels) - 1):
            cv2.imwrite( "midia/" + name + str( k ) + ".png", self.levels[k].getLevel( img, parameters ) )
        

#How to instantiate and use this class
#if __name__ == '__main__':
#    params = Parameters('params.yaml')
#    img = cv2.imread('../../box_in_scene.png')
#    model = cv2.imread('../../box.png')
#    rows,cols = img.shape[:2]
#    u = [ cols, rows ]
#    fovea = Fovea( u, params )
#    output = fovea.foveatedImage( img, params )
#    cv2.imshow( "foveated image", output )
#    cv2.waitKey( 0 )
#    fovea.foveatedFeatures( img, params )
#    fovea.features.matchingFeatures( model, params ) 
#    cv2.destroyAllWindows()
