#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file main.py

\brief This file executable to fovea structure

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
import sys
import numpy as np
from parameters import Parameters
from fovea import Fovea

class Structure:
    def __init__( self, parameters, fovea ):
        self.parameters = parameters
        self.fovea = fovea

    def update( self, position ):
        #print( str(position[0]) + ", " + str(position[1]) )
        self.parameters.setFovea( position )
        self.fovea.updateFovea( self.parameters )


def callback( event, x, y, flags, params ):
    '''
    \fn callback( event, x, y, flags, params )
    
    \brief This function is responsible for controlling all actions with the user's mouse
    '''
    #print( str(x) + ", " + str(y) )
    position = [ x, y ]
    params.update( position )


def main():
    arguments = sys.argv[1:]
    print ('python main.py ''image'' ''file <yaml>'' ')
    print ('Example: python main.py ~/box_in_scene.png fovea.yaml')
    img = cv2.imread( arguments[0] )
    params = Parameters( arguments[1] )
    rows, cols = img.shape[:2]
    u = [ cols, rows ]
    fovea = Fovea( u, params )

    structure = Structure( params, fovea )
    
    winname = 'image foveated'
    cv2.namedWindow( winname )
    while ( True ):
        # bind the callback function to window
        cv2.setMouseCallback( winname, callback, structure )
        color = (255, 255, 255) # B G R
        output = fovea.foveatedImage( img, params, color )
        cv2.imshow( winname, output ) 
        key = cv2.waitKey( 2 )
        if ( key == ord('q') ):
            break

        # Controls from fovea
        if ( key == ord('d') ):
            params.w[0] = min( params.u[0], params.w[0] + 10 )
        if ( key == ord('a') ):
            params.w[0] = max( 1, params.w[0] - 10 )
        if ( key == ord('c') ):
            params.w[1] = min( params.u[1] - 1, params.w[1] + 10 )
        if ( key == ord('z') ):
            params.w[1] = max( 1, params.w[1] - 10 )
        
        if ( key == ord('s') ):
            fovea.saveFoveatedImage( "output", output, params )

        if ( key == ord('l') ):
            fovea.saveLevels( "teste", img, params )

    cv2.destroyAllWindows()


if __name__ == '__main__':
    main()

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

