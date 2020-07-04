#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
\file main2.py

\brief This file executable to multifovea structure

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
from multifovea import Multifovea

class Structure:
    def __init__( self, parameters, multifovea ):
        self.parameters = parameters
        self.multifovea = multifovea
        self.indexFovea = 0

    def updateFovea( self ):
        if ( self.indexFovea + 1 <= int( len(self.parameters.f)/2) ):
            self.indexFovea += 1
        else:
            self.indexFovea = 0

    def update( self, position ):
        #print( str(position[0]) + ", " + str(position[1]) )
        self.parameters.f[2*self.indexFovea] = position[0]
        self.parameters.f[2*self.indexFovea+1] = position[1]
        self.parameters.fixFovea()
        self.multifovea.updateFovea( self.indexFovea, self.parameters )


def callback( event, x, y, flags, structure ):
    '''
    \fn callback( event, x, y, flags, params )
    
    \brief This function is responsible for controlling all actions with the user's mouse
    '''
    #print( str(x) + ", " + str(y) )
    position = [ x - (structure.parameters.u[0]/2), y - (structure.parameters.u[1]/2) ]
    structure.update( position )


def main():
    arguments = sys.argv[1:]
    if ( len( arguments ) != 3 ):
        print ('python main.py ''< model >'' ''< scene >'' ''<yaml>'' ')
        print ('Example: python main.py ~/box.png ~/box_in_scene.png fovea.yaml')
    else:
        model = cv2.imread( arguments[0] )
        scene = cv2.imread( arguments[1] )
        params = Parameters( arguments[2] )
        rows, cols = scene.shape[:2]
        u = [ cols, rows ]
        multifovea = Multifovea( u, params )
        
        structure = Structure( params, multifovea )
        
        winname = 'image foveated'
        cv2.namedWindow( winname )
        while ( True ):
            # bind the callback function to window
            cv2.setMouseCallback( winname, callback, structure )
            output = multifovea.multifoveatedImage( scene, params )
            cv2.imshow( winname, output ) 
            key = cv2.waitKey( 2 )
            if ( key == ord('q') ):
                break
            
            #'''
            # Controls from fovea
            if ( key == ord('d') ):
                params.w[0] = min( params.u[0], params.w[0] + 10 )
            if ( key == ord('a') ):
                params.w[0] = max( 1, params.w[0] - 10 )
            if ( key == ord('c') ):
                params.w[1] = min( params.u[1] - 1, params.w[1] + 10 )
            if ( key == ord('z') ):
                params.w[1] = max( 1, params.w[1] - 10 )
            if ( key == ord('m') ):
                structure.updateFovea()
            #'''
            
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

