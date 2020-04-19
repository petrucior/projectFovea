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
        p = parameters
        self.foveas = [] # cleaning foveas
        # Iteration 2 to 2 [(fx1, fy1), (fx2, fy2), ...]
        for k in range( 0, len(parameters.f), 2):
            px = [ parameters.f[k], parameters.f[k+1] ]
            p.setFovea( px )
            fovea = Fovea( u, p )
            self.foveas.append( fovea )
    
    def updateFovea( self, index, parameters ):
        '''
        \fn updateFovea( index, parameters )
        
        \brief This method update the fovea structure
        
        \param parameters - Parameters of fovea structure
        '''
        self.foveas[index].updateFovea( parameters )
    
    
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
