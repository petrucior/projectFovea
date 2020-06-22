#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file fovea.py

\brief This file contains the prototype of multifovea.

\author
Petrucio Ricardo Tavares de Medeiros \n
Universidade Federal do Rio Grande do Norte \n 
Departamento de Computacao e Automacao Industrial \n
petrucior at ufrn (dot) edu (dot) br

\version 0.1
\date January 2020

\note
This code demands the instalation of matplotlib 3d library:
pip install mplot3d-dragger ( see https://pypi.org/project/mplot3d-dragger/ )
 
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
from fovea import Fovea
import copy

class Multifovea :
    '''
    \class Multifovea
    
    \brief This class implements the Multifovea TAD
    '''
    
    def __init__( self, u, parameters ):
        '''
        \fn Multifovea( u, parameters )
        
        \brief Constructor default
        
        \param u - Size of image
        \param parameters - Parameters of fovea structure
        '''
        parameters.updateParameterSizeImage( u )
        p = copy.copy(parameters)
        self.foveas = [] # cleaning foveas
        # pixels
        px = []
        for i in range( 0, len(parameters.f), 2 ):
            pixel = [ parameters.f[i], parameters.f[i+1] ]
            px.append( pixel )
        # colors
        colors = []
        for c in range( 0, len(parameters.colors), 3):
            color = ( parameters.colors[c], parameters.colors[c+1], parameters.colors[c+2] )
            colors.append( color )
        # Iteration 2 to 2 [(fx1, fy1), (fx2, fy2), ...]
        for i in range( 0, int(len(parameters.f)/2) ):
            p.f = px[i]
            p.colors = colors[i]
            #p.updateParameterSizeImage( u )
            fovea = Fovea( u, p )
            fovea.updateFovea( p )
            self.foveas.append( fovea )
    
    def updateFovea( self, index, parameters ):
        '''
        \fn updateFovea( index, parameters )
        
        \brief This method update the fovea structure
        
        \param parameters - Parameters of fovea structure
        '''
        p = copy.copy(parameters)
        # pixels
        px = []
        for i in range( 0, len(parameters.f), 2 ):
            pixel = [ parameters.f[i], parameters.f[i+1] ]
            px.append( pixel )
        # colors
        colors = []
        for c in range( 0, len(parameters.colors), 3):
            color = ( parameters.colors[c], parameters.colors[c+1], parameters.colors[c+2] )
            colors.append( color )
        p.f = px[index]
        p.colors = colors[index]
        self.fovea[index].updateFovea( p )
        
        
    def multifoveatedFeatures( self, img, parameters ):
        '''
        \fn multifoveatedFeatures( parameters )
        
        \brief Compute keypoints and descriptors of levels
        
        \param img - Image to be foveated
        \param parameters - Parameters of fovea structure
        '''
        p = copy.copy(parameters)
        # pixels
        px = []
        for i in range( 0, len(parameters.f), 2 ):
            pixel = [ parameters.f[i], parameters.f[i+1] ]
            px.append( pixel )
        # colors
        colors = []
        for c in range( 0, len(parameters.colors), 3):
            color = ( parameters.colors[c], parameters.colors[c+1], parameters.colors[c+2] )
            colors.append( color )
        for f in range( 0, len(self.foveas) ):
            if ( parameters.typeFovea == 0 ): # MRMF
                p.f = px[f]
                p.colors = colors[f]
                self.foveas[f].foveatedFeatures( img, p )
            else: # MMF
                print("coming soon . . .")
            
            
    def multifoveatedMatching( self, model, parameters ):
        '''
        \fn foveatedMatching( model, parameters )
        
        \brief Maching between multiresolution domain and model 
        
        \param model - Template that represent the visual stimulli
        \param parameters - Parameters of fovea structure
        '''
        p = copy.copy(parameters)
        # pixels
        px = []
        for i in range( 0, len(parameters.f), 2 ):
            pixel = [ parameters.f[i], parameters.f[i+1] ]
            px.append( pixel )
        # colors
        colors = []
        for c in range( 0, len(parameters.colors), 3):
            color = ( parameters.colors[c], parameters.colors[c+1], parameters.colors[c+2] )
            colors.append( color )
        for f in range( 0, len(self.foveas) ):
            p.f = px[f]
            p.colors = colors[f]
            self.foveas[f].foveatedMatching( model, p )


    def multifoveatedImage( self, img, parameters ):
        '''
        \fn multifoveatedImage( img, parameters )
        
        \brief This function builds the multifocused image.
    
        \param img - Image to be foveated
        \param parameters - Parameters of fovea structure
        '''
        p = copy.copy(parameters)
        imgMultiFoveated = img.copy()
        # pixels
        px = []
        for i in range( 0, len(parameters.f), 2 ):
            pixel = [ parameters.f[i], parameters.f[i+1] ]
            px.append( pixel )
        # colors
        colors = []
        for c in range( 0, len(parameters.colors), 3):
            color = ( parameters.colors[c], parameters.colors[c+1], parameters.colors[c+2] )
            colors.append( color )
        for k in range( 0, parameters.m + 1 ):
            for f in range( 0, len(self.foveas) ):
                p.f = px[f]
                p.colors = colors[f]
                imgLevel = self.foveas[f].levels[k].getLevel( img, p )
                pi = self.foveas[f].mapLevel2Image( k, p, [0, 0] )
                pf = self.foveas[f].mapLevel2Image( k, p, p.w )
                if ( k < p.m ):
                    imgMultiFoveated[ int(pi[1]):int(pf[1]), int(pi[0]):int(pf[0]) ] = cv2.resize( imgLevel, ( int(pf[0]) - int(pi[0]), int(pf[1]) - int(pi[1]) ) )
                else:
                    imgMultiFoveated[ int(pi[1]):int(pf[1]), int(pi[0]):int(pf[0]) ] = imgLevel
                imgMultiFoveated = cv2.rectangle(imgMultiFoveated, (int(pi[0]), int(pi[1])), (int(pf[0]-1), int(pf[1]-1)), p.colors )
        return imgMultiFoveated
        
    
#How to instantiate and use this class
'''
if __name__ == '__main__':
    params = Parameters('params.yaml')
    model = cv2.imread('../../box.png')
    scene = cv2.imread('../../box_in_scene.png')
    rows,cols = scene.shape[:2]
    u = [ cols, rows ]
    multifovea = Multifovea( u, params )
    output = multifovea.multifoveatedImage( scene, params )
    cv2.imshow( "multifoveated image", output )
    cv2.waitKey( 0 )
    cv2.destroyAllWindows()
    multifovea.multifoveatedFeatures( scene, params )
    multifovea.multifoveatedMatching( model, params )
'''

